/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.cpp
 *
 * Description: This file contains the implementation of the class MiniStep.
 *****************************************************************************/

#include "MiniStep.h"
#include "model/variables/DependentVariable.h"
#include "model/ml/Chain.h"

namespace siena
{

/**
 * Constructs a new ministep.
 * @param[in] ego the actor making the change
 * @param[in] variableName the name of the dependent variable to be changed
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
MiniStep::MiniStep(int ego, string variableName, int difference)
{
	this->lego = ego;
	this->lvariableName = variableName;
	this->ldifference = difference;
	this->lpChain = 0;
	this->llogOptionSetProbability = 0;
	this->llogChoiceProbability = 0;
	this->lreciprocalRate = 0;
	this->lpPrevious = 0;
	this->lpNext = 0;
	this->lindex = -1;
	this->ldiagonalIndex = -1;
}


/**
 * Deallocates this ministep.
 */
MiniStep::~MiniStep()
{
}


/**
 * Returns the name of the dependent variable that this ministep is changing.
 */
string MiniStep::variableName() const
{
	return this->lvariableName;
}


/**
 * Stores the owner chain of this ministep.
 */
void MiniStep::pChain(Chain * pChain)
{
	this->lpChain = pChain;
}


/**
 * Returns the owner chain of this ministep.
 */
Chain * MiniStep::pChain() const
{
	return this->lpChain;
}


/**
 * Stores the log probability of choosing the option set of this ministep,
 * given the state just before this ministep.
 */
void MiniStep::logOptionSetProbability(double probability)
{
	this->llogOptionSetProbability = probability;
}


/**
 * Stores the log probability of making this ministep,
 * given that a ministep of the same option set will be made.
 */
void MiniStep::logChoiceProbability(double probability)
{
	this->llogChoiceProbability = probability;
}


/**
 * Stores the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
void MiniStep::reciprocalRate(double value)
{
	if (this->lpChain)
	{
		this->lpChain->onReciprocalRateChange(this, value);
	}

	this->lreciprocalRate = value;
}


/**
 * Stores the pointer to the previous ministep in the same chain.
 */
void MiniStep::pPrevious(MiniStep * pMiniStep)
{
	this->lpPrevious = pMiniStep;
}


/**
 * Stores the pointer to the next ministep in the same chain.
 */
void MiniStep::pNext(MiniStep * pMiniStep)
{
	this->lpNext = pMiniStep;
}


/**
 * Changes the given dependent variable according to this ministep.
 */
void MiniStep::makeChange(DependentVariable * pVariable)
{
	// Nothing in the base class.
}


/**
 * Returns if this ministep is diagonal, namely, it does not change
 * the dependent variables. Dummy ministeps are not considered diagonal.
 */
bool MiniStep::diagonal() const
{
	return false;
}

}
