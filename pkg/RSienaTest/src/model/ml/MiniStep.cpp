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
 * @param[in] variableId the ID of the dependent variable to be changed
 * @param[in] ego the actor making the change
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
MiniStep::MiniStep(int variableId, int ego, int difference)
{
	this->lego = ego;
	this->lvariableId = variableId;
	this->ldifference = difference;
	this->lpChain = 0;
	this->llogOptionSetProbability = 0;
	this->llogChoiceProbability = 0;
	this->lreciprocalRate = 0;
	this->lpPrevious = 0;
	this->lpNext = 0;
	this->lpPreviousWithSameOption = 0;
	this->lpNextWithSameOption = 0;
	this->lindex = -1;
	this->ldiagonalIndex = -1;
	this->lorderingKey = 0;
}


/**
 * Deallocates this ministep.
 */
MiniStep::~MiniStep()
{
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
 * Stores the pointer to the previous ministep having the same option as this
 * ministep.
 */
void MiniStep::pPreviousWithSameOption(MiniStep * pMiniStep)
{
	this->lpPreviousWithSameOption = pMiniStep;
}


/**
 * Stores the pointer to the next ministep having the same option as this
 * ministep.
 */
void MiniStep::pNextWithSameOption(MiniStep * pMiniStep)
{
	this->lpNextWithSameOption = pMiniStep;
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
