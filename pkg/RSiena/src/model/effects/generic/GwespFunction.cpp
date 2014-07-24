/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespFunction. Modified by Nynke Niezink, 10/02/14.
 *****************************************************************************/

#include "GwespFunction.h"
#include "network/Network.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"
#include <math.h>

namespace siena
{

/**
 * Constructor.
 */
GwespFunction::GwespFunction(string networkName,
	EgocentricConfigurationTable *	(NetworkCache::*pTable)() const,
	double parameter) :
	NetworkAlterFunction(networkName)
{
	this->lparameter = parameter;
	this->lweight = -0.01 * this->lparameter;
	this->lpTable = pTable;
}

/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void GwespFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpInitialisedTable = (*this->pNetworkCache().*lpTable)();

	// initialize the vector with weights for the GWESP statistic
	// this is done several times during one estimation run (not elegant,
	// but not computationally burdensome either)
	double pow = 1;
	int n = this->pNetwork()->n();
	this->lcumulativeWeight.resize(n); // default values 0
	for (int i = 1; i < n; i++)
	{
		pow *= (1 - exp(this->lweight));
		this->lcumulativeWeight[i] = exp(-(this->lweight)) * (1 - pow);
	}
}

/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double GwespFunction::value(int alter)
{
	int statistic = lpInitialisedTable->get(alter);
	return this->lcumulativeWeight[statistic];
}


}
