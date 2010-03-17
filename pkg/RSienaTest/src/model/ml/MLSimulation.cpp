#include <R.h>
#include "MLSimulation.h"
#include "utils/Random.h"
#include "model/SimulationActorSet.h"
#include "model/variables/DependentVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/ml/Chain.h"
#include "model/ml/MiniStep.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/BehaviorChange.h"

namespace siena
{

MLSimulation::MLSimulation(Data * pData, Model * pModel) :
	EpochSimulation(pData, pModel)
{
	this->lpChain = new Chain(pData);
}


MLSimulation::~MLSimulation()
{
	delete this->lpChain;
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void MLSimulation::connect(int period)
{
	this->lpChain->connect(period);
	this->initialize(period);
}


/**
 * Updates the probabilities for a range of ministeps of the given chain
 * (including the end-points of the range).
 */
void MLSimulation::updateProbabilities(Chain * pChain,
	MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	// Initialize the variables as of the beginning of the period
    this->initialize(pChain->period());

    // Apply the ministeps before the first ministep of the required range
    // to derive the correct state before that ministep.

    this->executeMiniSteps(pChain->pFirst()->pNext(), pFirstMiniStep);

    bool done = false;
    MiniStep * pMiniStep = pFirstMiniStep;

    while (!done)
    {
    	DependentVariable * pVariable =
    		this->lvariableMap[pMiniStep->variableName()];
    	this->calculateRates();
    	double rate = pVariable->rate(pMiniStep->ego());
    	double probability = pVariable->probability(pMiniStep);
    	double reciprocalTotalRate = 1 / this->totalRate();

    	pMiniStep->reciprocalRate(reciprocalTotalRate);
    	pMiniStep->logOptionSetProbability(log(rate * reciprocalTotalRate));
    	pMiniStep->logChoiceProbability(log(probability));
    	pMiniStep->makeChange(pVariable);

    	if (pMiniStep == pLastMiniStep)
    	{
    		done = true;
    	}
    	else
    	{
    		pMiniStep = pMiniStep->pNext();
    	}
    }
}


/**
 * Executes the given subsequence of ministeps excluding the last
 * ministep.
 */
void MLSimulation::executeMiniSteps(MiniStep * pFirstMiniStep,
	MiniStep * pLastMiniStep)
{
	MiniStep * pMiniStep = pFirstMiniStep;

	while (pMiniStep != pLastMiniStep)
	{
		DependentVariable * pVariable =
			this->lvariableMap[pMiniStep->variableName()];
		pMiniStep->makeChange(pVariable);

		pMiniStep = pMiniStep->pNext();
	}
}


void MLSimulation::setStateBefore(MiniStep * pMiniStep)
{
	this->resetVariables();
	this->executeMiniSteps(this->lpChain->pFirst()->pNext(), pMiniStep);
}


void MLSimulation::resetVariables()
{
	// Initialize each dependent variable

	for (unsigned i = 0; i < this->rVariables().size(); i++)
	{
		this->rVariables()[i]->initialize(this->period());
	}
}


// ----------------------------------------------------------------------------
// Section: Metropolis-Hastings steps
// ----------------------------------------------------------------------------

bool MLSimulation::insertDiagonalMiniStep()
{
	bool accept = false;
	MiniStep * pMiniStep = this->lpChain->randomMiniStep();
	this->setStateBefore(pMiniStep);
	this->calculateRates();
	DependentVariable * pVariable = this->chooseVariable();
	int i = this->chooseActor(pVariable);
	BehaviorVariable * pBehaviorVariable =
		dynamic_cast<BehaviorVariable *>(pVariable);

	if (pVariable->pActorSet()->active(i) &&
		(!pBehaviorVariable || !pBehaviorVariable->structural(i)))
	{
		MiniStep * pNewMiniStep = 0;

		if (pBehaviorVariable)
		{
			pNewMiniStep = new BehaviorChange(i, pVariable->name(), 0);
		}
		else
		{
			pNewMiniStep = new NetworkChange(i, i, pVariable->name(), 0);
		}

		double rr = 1 / this->totalRate();
		pMiniStep->reciprocalRate(rr);
		pMiniStep->logOptionSetProbability(log(pVariable->rate(i) * rr));
		pMiniStep->logChoiceProbability(
			log(pVariable->probability(pNewMiniStep)));

		double kappaFactor;

		if (this->lsimpleRates)
		{
			kappaFactor = 1 / (rr * (this->lpChain->ministepCount() + 1));
		}
		else
		{
			double sigma2 = this->lpChain->sigma2();
			double mu = this->lpChain->mu();

			kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu - rr) * (1 - mu - rr) / (2 * (sigma2 + rr * rr)));
		}

		this->lproposalProbability =
			kappaFactor * pVariable->probability(pNewMiniStep) *
				this->lpChain->ministepCount() *
				this->lcancelDiagonalProbability /
			((this->lpChain->diagonalMinistepCount() + 1) *
				this->linsertDiagonalProbability);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		if (nextDouble() < this->lproposalProbability)
		{
			accept = true;
			this->lpChain->insertBefore(pNewMiniStep, pMiniStep);
		}
	}

	return accept;
}


bool MLSimulation::cancelDiagonalMiniStep()
{
	bool accept = false;

	MiniStep * pMiniStep = this->lpChain->randomDiagonalMiniStep();
	double rr = pMiniStep->reciprocalRate();

	double kappaFactor;

	if (this->lsimpleRates)
	{
		kappaFactor = rr * this->lpChain->ministepCount();
	}
	else
	{
		double sigma2 = this->lpChain->sigma2();
		double mu = this->lpChain->mu();

		kappaFactor = sqrt(sigma2 / (sigma2 + rr * rr)) *
			exp((1 - mu) * (1 - mu) / (2 * sigma2) -
				(1 - mu + rr) * (1 - mu + rr) / (2 * (sigma2 - rr * rr)));
	}

	this->lproposalProbability =
		kappaFactor * exp(-pMiniStep->logChoiceProbability()) *
			this->lpChain->diagonalMinistepCount() *
			this->linsertDiagonalProbability /
		((this->lpChain->ministepCount() + 1) *
			this->lcancelDiagonalProbability);

	if (this->lproposalProbability > 1)
	{
		this->lproposalProbability = 1;
	}

	if (nextDouble() < this->lproposalProbability)
	{
		accept = true;
		this->lpChain->remove(pMiniStep);
	}

	return accept;
}


bool MLSimulation::permute(int c0)
{
	if (this->lpChain->ministepCount() <= 2)
	{
		return false;
	}

	bool accept = false;

	MiniStep * pMiniStepA = this->lpChain->randomMiniStep();

	while (pMiniStepA == this->lpChain->pLast())
	{
		pMiniStepA = this->lpChain->randomMiniStep();
	}

	vector<MiniStep *> interval;
	MiniStep * pMiniStep = pMiniStepA;

	while ((int) interval.size() < c0 && pMiniStep != this->lpChain->pLast())
	{
		interval.push_back(pMiniStep);
		pMiniStep = pMiniStep->pNext();
	}

	if (interval.size() <= 1)
	{
		return false;
	}

	MiniStep * pNextMiniStep = pMiniStep;

	permuteVector(interval);
	this->setStateBefore(pMiniStepA);
	bool valid = true;
	double sumlprob = 0;
	double sumlprob_new = 0;
	double mu_new = this->lpChain->mu();
	double sigma2_new = this->lpChain->sigma2();
	double * newReciprocalRate = new double[interval.size()];
	double * newOptionSetProbability = new double[interval.size()];
	double * newChoiceProbability = new double[interval.size()];

	for (unsigned i = 0; i < interval.size() && valid; i++)
	{
		pMiniStep = interval[i];
		DependentVariable * pVariable =
			this->lvariableMap[pMiniStep->variableName()];

		if (!pVariable->validMiniStep(pMiniStep))
		{
			valid = false;
		}
		else
		{
			sumlprob += pMiniStep->logChoiceProbability() +
				pMiniStep->logOptionSetProbability();
			double rrOld = pMiniStep->reciprocalRate();

			if (!this->simpleRates())
			{
				mu_new -= rrOld;
				sigma2_new -= rrOld * rrOld;
			}

			this->calculateRates();
			double rr = 1 / this->totalRate();
			double lospr =
				log(pVariable->rate(pMiniStep->ego()) *
					newReciprocalRate[i]);
			double lcpr = log(pVariable->probability(pMiniStep));

			sumlprob_new += lospr + lcpr;

			if (!this->simpleRates())
			{
				mu_new += rr;
				sigma2_new += rr * rr;
			}

			pMiniStep->makeChange(pVariable);

			newReciprocalRate[i] = rr;
			newOptionSetProbability[i] = lospr;
			newChoiceProbability[i] = lcpr;
		}
	}

	if (valid)
	{
		double kappaFactor = 1;

		if (!this->simpleRates())
		{
			double sigma2 = this->lpChain->sigma2();
			double mu = this->lpChain->mu();

			kappaFactor = sqrt(sigma2 / sigma2_new) *
				exp((1 - mu) * (1 - mu) / (2 * sigma2) -
					(1 - mu_new) * (1 - mu_new) / (2 * sigma2_new));
		}

		this->lproposalProbability =
			kappaFactor * exp(sumlprob_new - sumlprob);

		if (this->lproposalProbability > 1)
		{
			this->lproposalProbability = 1;
		}

		if (nextDouble() < this->lproposalProbability)
		{
			accept = true;

			for (unsigned i = 0; i < interval.size(); i++)
			{
				pMiniStep = interval[i];

				this->lpChain->remove(pMiniStep);

				pMiniStep->reciprocalRate(newReciprocalRate[i]);
				pMiniStep->logOptionSetProbability(
					newOptionSetProbability[i]);
				pMiniStep->logChoiceProbability(newChoiceProbability[i]);
			}

			for (unsigned i = 0; i < interval.size(); i++)
			{
				this->lpChain->insertBefore(interval[i], pNextMiniStep);
			}
		}
	}

	delete[] newReciprocalRate;
	delete[] newOptionSetProbability;
	delete[] newChoiceProbability;

	return accept;
}


// ----------------------------------------------------------------------------
// Section: Accessors
// ----------------------------------------------------------------------------

/**
 * Returns the proposal probability of the last Metropolis-Hastings step.
 */
double MLSimulation::proposalProbability() const
{
	return this->lproposalProbability;
}


/**
 * Stores if simple rates should be used in simulations.
 */
void MLSimulation::simpleRates(bool flag)
{
	this->lsimpleRates = flag;
}


/**
 * Returns if simple rates should be used in simulations.
 */
bool MLSimulation::simpleRates() const
{
	return this->lsimpleRates;
}


/**
 * Stores the probability associated with the insertDiagonalMiniStep
 * operation.
 */
void MLSimulation::insertDiagonalProbability(double probability)
{
	this->linsertDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the insertDiagonalMiniStep
 * operation.
 */
double MLSimulation::insertDiagonalProbability() const
{
	return this->linsertDiagonalProbability;
}


/**
 * Stores the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
void MLSimulation::cancelDiagonalProbability(double probability)
{
	this->lcancelDiagonalProbability = probability;
}


/**
 * Returns the probability associated with the cancelDiagonalMiniStep
 * operation.
 */
double MLSimulation::cancelDiagonalProbability() const
{
	return this->lcancelDiagonalProbability;
}

}
