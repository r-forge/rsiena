/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: SameCovariateInTiesFunction.h
 *
 * Description: This file contains the definition of the
 * SameCovariateTwoPathFunction class.
 *****************************************************************************/

#ifndef SAMECOVARIATEINTIESFUNCTION_H_
#define SAMECOVARIATEINTIESFUNCTION_H_

#include "CovariateNetworkAlterFunction.h"
#include "model/effects/NetworkEffect.h"

namespace siena
{

class SameCovariateInTiesFunction: public CovariateNetworkAlterFunction
{
public:
	SameCovariateInTiesFunction(string networkName,
		string covariateName, bool excludeMissing);

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

		virtual double value(int alter);

private:
	bool lexcludeMissing;
};

}

#endif /* SAMECOVARIATEINTIESFUNCTION_H_ */
