/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: BehaviorVariable.cpp
 *
 * Description: This file contains the implementation of the
 * BehaviorVariable class.
 *****************************************************************************/
#include <cstdlib>
#include <cmath>
#include <string>
#include <stdexcept>
#include "data/ActorSet.h"
#include "utils/Random.h"
#include "BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"
#include "model/EpochSimulation.h"
#include "model/tables/Cache.h"
#include "model/variables/DependentVariable.h"
#include "model/Model.h"
#include "model/effects/BehaviorEffect.h"
#include "model/EffectInfo.h"
#include "model/SimulationActorSet.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"

namespace siena
{

/**
 * Creates a new behavior variable for the given observed data.
 * @param pModel the owner model of this variable
 */
BehaviorVariable::BehaviorVariable(BehaviorLongitudinalData * pData,
	EpochSimulation * pSimulation) :
		DependentVariable(pData->name(),
			pData->pActorSet(),
			pSimulation)
{
	this->lpData = pData;
	this->lvalues = new int[this->n()];
	this->levaluationEffectContribution = new double * [3];
	this->lendowmentEffectContribution = new double * [3];
	this->lprobabilities = new double[3];

	for (int i = 0; i < 3; i++)
	{
		this->levaluationEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
		this->lendowmentEffectContribution[i] =
			new double[pSimulation->pModel()->rEvaluationEffects(pData->name()).size()];
	}
}


/**
 * Deallocates this variable object.
 */
BehaviorVariable::~BehaviorVariable()
{
	delete[] this->lvalues;

	this->lpData = 0;
	this->lvalues = 0;
	delete[] this->lprobabilities;
	// Delete arrays of contributions

	for (int i = 0; i < 3; i++)
	{
		delete[] this->levaluationEffectContribution[i];
		delete[] this->lendowmentEffectContribution[i];
	}

	delete[] this->levaluationEffectContribution;
	delete[] this->lendowmentEffectContribution;

	this->levaluationEffectContribution = 0;
	this->lendowmentEffectContribution = 0;
	this->lprobabilities = 0;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the second dimension of this variable, namely, how many values
 * correspond to each actor. This number is 1 for behavior variables.
 */
int BehaviorVariable::m() const
{
	return 1;
}


/**
 * Returns the longitudinal data object this variable is based on.
 */
LongitudinalData * BehaviorVariable::pData() const
{
	return this->lpData;
}


/**
 * Returns if the observed start value of the given actor is missing
 * for the current period.
 */
bool BehaviorVariable::missingStartValue(int actor) const
{
	return this->lpData->missing(this->period(), actor);
}


/**
 * Returns the current value on this behavior for the given actor.
 */
int BehaviorVariable::value(int actor) const
{
	return this->lvalues[actor];
}


/**
 * Stores the current value on this behavior for the given actor.
 */
void BehaviorVariable::value(int actor, int newValue)
{
	this->lvalues[actor] = newValue;
}


/**
 * Returns the current value on this behavior for the given actor, which is
 * centered around the overall mean of the observed values.
 */
double BehaviorVariable::centeredValue(int actor) const
{
	return this->lvalues[actor] - this->lpData->overallMean();
}


/**
 * Returns if the behavior is structurally determined for the given actor
 * at the current period.
 */
bool BehaviorVariable::structural(int actor) const
{
	return this->lpData->structural(this->period(), actor);
}


/**
 * Returns the centered similarity of the given actors.
 */
double BehaviorVariable::similarity(int i, int j) const
{
	return this->lpData->similarity(this->lvalues[i], this->lvalues[j]);
}


/**
 * Returns the array of current values for this behavior variable.
 */
const int * BehaviorVariable::values() const
{
	return this->lvalues;
}


/**
 * Returns the range of observed values for this behavior variable.
 */
int BehaviorVariable::range() const
{
	return this->lpData->range();
}

/**
 * Returns the similarity mean for this behavior variable.
 */
double BehaviorVariable::similarityMean() const
{
	return this->lpData->similarityMean();
}

// ----------------------------------------------------------------------------
// Section: Initialization at the beginning of a period
// ----------------------------------------------------------------------------

/**
 * Initializes this variable as of the beginning of the given period.
 */
void BehaviorVariable::initialize(int period)
{
	DependentVariable::initialize(period);

	// Copy the values from the corresponding observation.

	for (int i = 0; i < this->n(); i++)
	{
		this->lvalues[i] = this->lpData->value(period, i);
	}
}


// ----------------------------------------------------------------------------
// Section: Composition change
// ----------------------------------------------------------------------------

/**
 * Sets leavers values back to the value at the start of the simulation.
 */
void BehaviorVariable::setLeaverBack(const SimulationActorSet * pActorSet,
	int actor)
{
	if (pActorSet == this->pActorSet())
	{
		// Reset ties from the given actor to values at start

		for (int i = 0; i < this->n(); i++)
		{
			this->lvalues[actor] =	this->lpData->value(this->period(), actor);
		}
	}
}


// ----------------------------------------------------------------------------
// Section: Changing the behavior variable
// ----------------------------------------------------------------------------

/**
 * Simulates a change of the behavior according to the choice of the given
 * actor.
 */
void BehaviorVariable::makeChange(int actor)
{
	this->calculateProbabilities(actor);

	// Choose the change
	int difference = nextIntWithProbabilities(3, this->lprobabilities) - 1;

	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(difference + 1,
			this->lupPossible,
			this->ldownPossible);
	}

	if (this->pSimulation()->pModel()->needChain())
	{
		// insert ministep in chain
		BehaviorChange * pMiniStep =
			new BehaviorChange(this->lpData, actor, difference);
		this->pSimulation()->pChain()->insertBefore(pMiniStep,
			this->pSimulation()->pChain()->pLast());
		pMiniStep->logChoiceProbability(log(this->lprobabilities[difference
					+ 1]));
	}
	// Make the change

	if (difference != 0)
	{
		int oldValue = this->lvalues[actor];

		// Make the change
		this->lvalues[actor] += difference;

		// Update the distance from the observed data at the beginning of the
		// period. Actors with missing values at any of the endpoints of the
		// period don't contribute to the distance

		if (!this->lpData->missing(this->period(), actor) &&
			!this->lpData->missing(this->period() + 1, actor))
		{
			int observedValue = this->lpData->value(this->period(), actor);
			this->simulatedDistance(this->simulatedDistance() +
				std::abs(this->lvalues[actor] - observedValue) -
				std::abs(oldValue - observedValue));
		}
	}
}


/**
 * Calculates the probabilities of each possible change.
 */
void BehaviorVariable::calculateProbabilities(int actor)
{
	this->lupPossible = true;
	this->ldownPossible = true;

	// Calculate the probability for downward change

	if (this->lvalues[actor] > this->lpData->min() &&
		!this->lpData->upOnly(this->period()))
	{
		this->lprobabilities[0] =
			exp(this->totalEvaluationContribution(actor, -1) +
				this->totalEndowmentContribution(actor, -1));
	}
	else
	{
		this->lprobabilities[0] = 0;
		this->ldownPossible = false;
	}

	// No change means zero contribution, but exp(0) = 1
	this->lprobabilities[1] = 1;

	// Calculate the probability for upward change

	if (this->lvalues[actor] < this->lpData->max() &&
		!this->lpData->downOnly(this->period()))
	{
		this->lprobabilities[2] =
			exp(this->totalEvaluationContribution(actor, 1));
	}
	else
	{
		this->lprobabilities[2] = 0;
		this->lupPossible = false;
	}

	double sum = this->lprobabilities[0] +
		this->lprobabilities[1] +
		this->lprobabilities[2];
	this->lprobabilities[0] /= sum;
	this->lprobabilities[1] /= sum;
	this->lprobabilities[2] /= sum;
}


/**
 * Returns the total contribution of all effects in the given function if
 * the behavior of the given actor is changed by the given amount.
 */
double BehaviorVariable::totalEvaluationContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEvaluationFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect =
			(BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		this->levaluationEffectContribution[difference+1][i] =
			thisContribution;
		contribution += pEffect->parameter() * thisContribution;
	}
	return contribution;
}

double BehaviorVariable::totalEndowmentContribution(int actor,
	int difference) const
{
	double contribution = 0;
	const Function * pFunction = this->pEndowmentFunction();

	for (unsigned i = 0; i < pFunction->rEffects().size(); i++)
	{
		BehaviorEffect * pEffect =
			(BehaviorEffect *) pFunction->rEffects()[i];
		double thisContribution =
			pEffect->calculateChangeContribution(actor, difference);
		this->lendowmentEffectContribution[difference+1][i] =
			thisContribution;
		contribution += pEffect->parameter() * thisContribution;
	}

	return contribution;
}
/**
 * Updates the scores for evaluation and endowment function effects according
 * to the current step in the simulation.
 */
void BehaviorVariable::accumulateScores(int difference,
	bool upPossible, bool downPossible) const
{
	for (unsigned i = 0;
		i < this->pEvaluationFunction()->rEffects().size();
		i++)
	{
		if (difference == 1) // no change, but not initialised
		{
			this->levaluationEffectContribution[difference][i] = 0;
		}
		Effect * pEffect = this->pEvaluationFunction()->rEffects()[i];
		double score = this->levaluationEffectContribution[difference][i];

		if (upPossible)
		{
			score -=
				this->levaluationEffectContribution[2][i] *
				this->lprobabilities[2];
		}

		if (downPossible)
		{
			score -=
				this->levaluationEffectContribution[0][i] *
				this->lprobabilities[0];
		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}

	for (unsigned i = 0;
		i < this->pEndowmentFunction()->rEffects().size();
		i++)
	{
		if (difference == 1) // no change, but not initialised
		{
			this->lendowmentEffectContribution[difference][i] = 0;
		}
		if (difference == 2) // up has no effect on endowment function
		{
			this->lendowmentEffectContribution[difference][i] = 0;
		}
		Effect * pEffect = this->pEndowmentFunction()->rEffects()[i];
		double score = this->lendowmentEffectContribution[difference][i];

		if (downPossible)
		{
			score -=
				this->lendowmentEffectContribution[0][i] *
					this->lprobabilities[0];

		}

		this->pSimulation()->score(pEffect->pEffectInfo(),
			this->pSimulation()->score(pEffect->pEffectInfo()) + score);
	}
}


// ----------------------------------------------------------------------------
// Section: Maximum likelihood related methods
// ----------------------------------------------------------------------------

/**
 * Calculates the probability of the given ministep assuming that the
 * ego of the ministep will change this variable.
 */
double BehaviorVariable::probability(MiniStep * pMiniStep)
{
	// Initialize the cache object for the current ego
	this->pSimulation()->pCache()->initialize(pMiniStep->ego());

	BehaviorChange * pBehaviorChange =
		dynamic_cast<BehaviorChange *>(pMiniStep);

	if (pBehaviorChange->difference() < -1 ||
		pBehaviorChange->difference() > 1)
	{
		throw invalid_argument("MiniStep difference out of range [-1,1].");
	}

	this->calculateProbabilities(pMiniStep->ego());
	if (this->pSimulation()->pModel()->needScores())
	{
		this->accumulateScores(pBehaviorChange->difference() + 1,
			this->lupPossible,
			this->ldownPossible);
	}
	if (this->pSimulation()->pModel()->needDerivatives())
	{
		this->accumulateDerivatives();
	}
	return this->lprobabilities[pBehaviorChange->difference() + 1];
}

/**
 * Updates the derivatives for evaluation and endowment function effects
 * according to the current miniStep in the chain.
 */
void BehaviorVariable::accumulateDerivatives() const
{
	int totalEvaluationEffects = this->pEvaluationFunction()->rEffects().size();
	int totalEndowmentEffects = this->pEndowmentFunction()->rEffects().size();
	int totalEffects = totalEvaluationEffects + totalEndowmentEffects;
	Effect * pEffect1;
	Effect * pEffect2;
	double derivative;
	double * product = new double[totalEffects];
	double contribution1 = 0.0;
	double contribution2 = 0.0;

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		product[effect1] = 0.0;

		if (effect1 < totalEvaluationEffects)
		{
			pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
		}
		else
		{
			pEffect1 = this->pEndowmentFunction()->rEffects()[effect1];
		}
		if (this->lupPossible)
		{
			if (effect1 < totalEvaluationEffects)
			{
				product[effect1] +=
					this->levaluationEffectContribution[2][effect1] *
					this->lprobabilities[2];
			}
			else
			{
				product[effect1] +=
					this->lendowmentEffectContribution[2][effect1] *
					this->lprobabilities[2];
			}
			//	Rprintf("%d %d %f\n", alter, effect1, product[effect1]);
		}
		if (this->ldownPossible)
		{
			if (effect1 < totalEvaluationEffects)
			{
				product[effect1] +=
					this->levaluationEffectContribution[0][effect1] *
					this->lprobabilities[0];
			}
			else
			{
				product[effect1] +=
					this->lendowmentEffectContribution[0][effect1] *
					this->lprobabilities[0];
			}
			//	Rprintf("%d %d %f\n", alter, effect1, product[effect1]);
		}
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			derivative = 0.0;
			if (effect2 <= totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[effect2];
			}

			if (this->lupPossible)
			{
				if (effect1 < totalEvaluationEffects)
				{
					contribution1 =
						this->levaluationEffectContribution[2][effect1];
				}
				else
				{
					contribution1 =
						this->lendowmentEffectContribution[2][effect1];
				}
				if (effect2 < totalEvaluationEffects)
				{
					contribution2 =
						this->levaluationEffectContribution[2][effect2];
				}
				else
				{
					contribution2 =
						this->lendowmentEffectContribution[2][effect2];
				}

				derivative -=
					contribution1 * contribution2 *	this->lprobabilities[2];
				//	Rprintf("deriv 2 %d %d %d %f %f %f %f\n", alter, effect1, effect2,
				//		derivative,
						//		this->levaluationEffectContribution[alter][effect1],
				//		this->levaluationEffectContribution[alter][effect2],
				//		this->lprobabilities[alter]);
			}

			if (this->ldownPossible)
			{
				if (effect1 < totalEvaluationEffects &&
					effect2 < totalEvaluationEffects)
				{
					contribution1 =
						this->levaluationEffectContribution[0][effect1];

					contribution2 =
						this->levaluationEffectContribution[0][effect2];
					derivative -=
						contribution1 * contribution2 *	this->lprobabilities[0];
				}
				//	Rprintf("deriv 2 %d %d %d %f %f %f %f\n", alter, effect1, effect2,
				//		derivative,
						//		this->levaluationEffectContribution[alter][effect1],
				//		this->levaluationEffectContribution[alter][effect2],
				//		this->lprobabilities[alter]);
			}
			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +	derivative);
		}
	}

	for (int effect1 = 0; effect1 < totalEffects; effect1++)
	{
		for (int effect2 = effect1; effect2 < totalEffects; effect2++)
		{
			if (effect1 < totalEvaluationEffects)
			{
				pEffect1 = this->pEvaluationFunction()->rEffects()[effect1];
			}
			else
			{
				pEffect1 = this->pEndowmentFunction()->rEffects()[effect1];
			}
			if (effect2 <= totalEvaluationEffects)
			{
				pEffect2 = this->pEvaluationFunction()->rEffects()[effect2];
			}
			else
			{
				pEffect2 = this->pEndowmentFunction()->rEffects()[effect2];
			}

			this->pSimulation()->derivative(pEffect1->pEffectInfo(),
				pEffect2->pEffectInfo(),
				this->pSimulation()->derivative(pEffect1->pEffectInfo(),
					pEffect2->pEffectInfo()) +
				product[effect1] * product[effect2]);
		}
	}
	delete[] product;

}

//	Rprintf("deriv %f\n", derivative;



/**
 * Returns whether applying the given ministep on the current state of this
 * variable would be valid with respect to all constraints.
 */
bool BehaviorVariable::validMiniStep(const MiniStep * pMiniStep) const
{
	bool valid = DependentVariable::validMiniStep(pMiniStep);

	if (valid && !pMiniStep->diagonal())
	{
		const BehaviorChange * pBehaviorChange =
			dynamic_cast<const BehaviorChange *>(pMiniStep);
		int i = pMiniStep->ego();
		int d = pBehaviorChange->difference();
		int newValue = this->lvalues[i] + d;

		if (newValue < this->lpData->min() || newValue > this->lpData->max())
		{
			valid = false;
		}
		else if (d > 0 && this->lpData->downOnly(this->period()))
		{
			valid = false;
		}
		else if (d < 0 && this->lpData->upOnly(this->period()))
		{
			valid = false;
		}
		else
		{
			valid = !this->lpData->structural(this->period(), i);
		}
	}

	return valid;
}


/**
 * Generates a random ministep for the given ego.
 */
MiniStep * BehaviorVariable::randomMiniStep(int ego)
{
	this->pSimulation()->pCache()->initialize(ego);
	this->calculateProbabilities(ego);
	int difference = nextIntWithProbabilities(3, this->lprobabilities) - 1;
	BehaviorChange * pMiniStep =
		new BehaviorChange(this->lpData, ego, difference);
	pMiniStep->logChoiceProbability(log(this->lprobabilities[difference + 1]));
	return pMiniStep;
}


/**
 * Returns if the observed value for the option of the given ministep
 * is missing at either end of the period.
 */
bool BehaviorVariable::missing(const MiniStep * pMiniStep) const
{
	return this->lpData->missing(this->period(), pMiniStep->ego()) ||
		this->lpData->missing(this->period() + 1, pMiniStep->ego());
}

/**
 * Returns if the given ministep is structurally determined in the period.
 */
bool BehaviorVariable::structural(const MiniStep * pMiniStep) const
{
	return this->lpData->structural(this->period(), pMiniStep->ego());
}

// ----------------------------------------------------------------------------
// Section: Properties
// ----------------------------------------------------------------------------

/**
 * Returns if this is a behavior variable.
 */
bool BehaviorVariable::behaviorVariable() const
{
	return true;
}

}
