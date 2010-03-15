/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MiniStep.h
 *
 * Description: This file contains the definition of the MiniStep class.
 *****************************************************************************/


#ifndef MINISTEP_H_
#define MINISTEP_H_

#include <string>

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class DependentVariable;
class Chain;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Defines a single ministep as part of a Chain object in the
 * Maximum-Likelihood calculations.
 */
class MiniStep
{
	friend class Chain;

public:
	MiniStep(int ego, string variableName, int difference);
	virtual ~MiniStep();

	inline int ego() const;
	string variableName() const;
	inline int difference() const;

	Chain * pChain() const;

	inline double logOptionSetProbability() const;
	void logOptionSetProbability(double probability);

	inline double logChoiceProbability() const;
	void logChoiceProbability(double probability);

	inline double reciprocalRate() const;
	void reciprocalRate(double value);

	inline MiniStep * pPrevious() const;
	inline MiniStep * pNext() const;
	void pPrevious(MiniStep * pMiniStep);
	void pNext(MiniStep * pMiniStep);

	inline int index() const;
	inline void index(int index);

	inline int diagonalIndex() const;
	inline void diagonalIndex(int index);

	virtual void makeChange(DependentVariable * pVariable);
	virtual bool diagonal() const;

private:
	void pChain(Chain * pChain);

	// The actor making the change
	int lego;

	// The name of the dependent variable to be changed
	string lvariableName;

	// The amount of change (+1,0,-1 for dichotomous variables)
	int ldifference;

	// The owner chain (0, if the ministep is not part of a chain)
	Chain * lpChain;

	// Log probability of choosing the option set of this ministep,
	// given the state just before this ministep.

	double llogOptionSetProbability;

	// Log probability of making this ministep, given that a ministep
	// of the same option set will be made.

	double llogChoiceProbability;

	// Reciprocal of aggregate (summed) rate function immediately
	// before this ministep.

	double lreciprocalRate;

	// Points to the previous ministep in the same chain
	MiniStep * lpPrevious;

	// Points to the next ministep in the same chain
	MiniStep * lpNext;

	// The index in the vector of ministeps of the owner chain.
	int lindex;

	// The index in the vector of diagonal ministeps of the owner chain.
	int ldiagonalIndex;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the ego of this ministep.
 */
int MiniStep::ego() const
{
	return this->lego;
}


/**
 * Returns the amount of change in this ministep
 * (+1,0,-1 for dichotomous variables).
 */
int MiniStep::difference() const
{
	return this->ldifference;
}


/**
 * Returns the log probability of choosing the option set of this ministep,
 * given the state just before this ministep.
 */
double MiniStep::logOptionSetProbability() const
{
	return this->llogOptionSetProbability;
}


/**
 * Returns the log probability of making this ministep,
 * given that a ministep of the same option set will be made.
 */
double MiniStep::logChoiceProbability() const
{
	return this->llogChoiceProbability;
}


/**
 * Returns the reciprocal of aggregate rate function immediately before
 * this ministep.
 */
double MiniStep::reciprocalRate() const
{
	return this->lreciprocalRate;
}


/**
 * Returns the previous ministep in the same chain.
 */
MiniStep * MiniStep::pPrevious() const
{
	return this->lpPrevious;
}


/**
 * Returns the next ministep in the same chain.
 */
MiniStep * MiniStep::pNext() const
{
	return this->lpNext;
}


/**
 * Returns the index of this ministep in the vector of ministeps of the
 * owner chain.
 */
int MiniStep::index() const
{
	return this->lindex;
}


/**
 * Stores the index of this ministep in the vector of ministeps of the
 * owner chain.
 */
void MiniStep::index(int index)
{
	this->lindex = index;
}


/**
 * Returns the index of this ministep in the vector of diagonal
 * ministeps of the owner chain (-1, if the ministep is not part of a
 * chain or is not diagonal).
 */
int MiniStep::diagonalIndex() const
{
	return this->ldiagonalIndex;
}


/**
 * Stores the index of this ministep in the vector of diagonal
 * ministeps of the owner chain (-1, if the ministep is not part of a
 * chain or is not diagonal).
 */
void MiniStep::diagonalIndex(int index)
{
	this->ldiagonalIndex = index;
}

}

#endif /* MINISTEP_H_ */
