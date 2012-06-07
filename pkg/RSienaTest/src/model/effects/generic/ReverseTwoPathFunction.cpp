/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: ReverseTwoPathFunction.cpp
 *
 * Description: This file contains the implementation of the class
 * ReverseTwoPathFunction.
 *****************************************************************************/

#include "ReverseTwoPathFunction.h"
#include "model/tables/NetworkCache.h"
#include "model/tables/EgocentricConfigurationTable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] networkName the name of the network variable this function is
 * associated with
 */
ReverseTwoPathFunction::ReverseTwoPathFunction(string networkName) :
	OneModeNetworkAlterFunction(networkName)
{
	this->lpTable = 0;
}


/**
 * Initializes this function.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void ReverseTwoPathFunction::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	OneModeNetworkAlterFunction::initialize(pData, pState, period, pCache);
	this->lpTable = this->pNetworkCache()->pReverseTwoPathTable();
}


/**
 * Returns the value of this function for the given alter. It is assumed
 * that the function has been initialized before and pre-processed with
 * respect to a certain ego.
 */
double ReverseTwoPathFunction::value(int alter)
{
	return this->lpTable->get(alter);
}


/**
 * Returns the value of this function as an integer.
 */
int ReverseTwoPathFunction::intValue(int alter)
{
	return this->lpTable->get(alter);
}

}
