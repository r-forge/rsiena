/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.h
 *
 * Description: This file contains the definition of the Chain class.
 *****************************************************************************/

#ifndef CHAIN_H_
#define CHAIN_H_

#include <vector>
#include <map>
#include "model/ml/Option.h"

using namespace std;

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class MiniStep;
class Data;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------


/**
 * Defines a sequence of ministeps, which is the basic structure for
 * Maximum-Likelihood calculations.
 */
class Chain
{
public:
	Chain(Data * pData);
	virtual ~Chain();

	// Chain modifications

	void clear();
	void insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep);
	void remove(MiniStep * pMiniStep);
	void connect(int period);

	void onReciprocalRateChange(const MiniStep * pMiniStep, double newValue);

	// Accessors

	int period() const;
	MiniStep * pFirst() const;
	MiniStep * pLast() const;
	int ministepCount() const;
	int diagonalMinistepCount() const;
	double mu() const;
	double sigma2() const;

	// Intervals

	int intervalLength(const MiniStep * pFirstMiniStep,
		const MiniStep * pLastMiniStep) const;

	// Same option related

	MiniStep * nextMiniStepForOption(Option & rOption,
		const MiniStep * pFirstMiniStep) const;

	// Random draws

	MiniStep * randomMiniStep() const;
	MiniStep * randomDiagonalMiniStep() const;
	MiniStep * randomMiniStep(MiniStep * pFirstMiniStep,
		MiniStep * pLastMiniStep) const;

private:
	void resetOrderingKeys();

	// A dummy first ministep in the chain
	MiniStep * lpFirst;

	// A dummy last ministep in the chain
	MiniStep * lpLast;

	// The underlying observed data
	Data * lpData;

	// The period of changes represented by this chain
	int lperiod;

	// Stores the ministeps in no particular order.
	// The first (dummy) ministep is not stored in this vector.

	vector<MiniStep *> lminiSteps;

	// Stores the diagonal ministeps in no particular order.
	vector<MiniStep *> ldiagonalMiniSteps;

	// Sum of reciprocal rates over all non-dummy ministeps.
	double lmu;

	// Sum of squared reciprocal rates over all non-dummy ministeps.
	double lsigma2;

	// Maps each option to its first ministep in the chain (if any)
	map<Option, MiniStep *> lfirstMiniStepPerOption;
};

}

#endif /* CHAIN_H_ */
