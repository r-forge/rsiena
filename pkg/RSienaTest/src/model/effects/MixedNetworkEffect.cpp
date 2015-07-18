/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: MixedNetworkEffect.cpp
 *
 * Description: This file contains the implementation of the
 * MixedNetworkEffect class.
 *****************************************************************************/

//#include <stdexcept>
//#include <R_ext/Print.h>

#include "MixedNetworkEffect.h"
#include "network/Network.h"
//#include "network/IncidentTieIterator.h"
//#include "data/NetworkLongitudinalData.h"
#include "model/State.h"
//#include "model/EpochSimulation.h"
#include "model/EffectInfo.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/Cache.h"
#include "model/tables/NetworkCache.h"
//#include "model/tables/EgocentricConfigurationTable.h"

using namespace std;

namespace siena
{

/**
 * Constructor.
 */
MixedNetworkEffect::MixedNetworkEffect(const EffectInfo * pEffectInfo,
	   string firstNetworkName,	string secondNetworkName) :
	NetworkEffect(pEffectInfo)
{
	this->lname1 = firstNetworkName;
	this->lname2 = secondNetworkName;
	this->lpFirstNetwork = 0;
	this->lpSecondNetwork = 0;
	this->lpTwoNetworkCache = 0;
	this->lpFirstNetworkCache = 0;
}

/**
 * Initializes this effect.
 * @param[in] pData the observed data
 * @param[in] pState the current state of the dependent variables
 * @param[in] period the period of interest
 * @param[in] pCache the cache object to be used to speed up calculations
 */
void MixedNetworkEffect::initialize(const Data * pData,
	State * pState,
	int period,
	Cache * pCache)
{
	NetworkEffect::initialize(pData, pState, period, pCache);
	this->lpFirstNetwork = pState->pNetwork(this->lname1);
	this->lpSecondNetwork = pState->pNetwork(this->lname2);
	this->lpTwoNetworkCache = pCache->pTwoNetworkCache(this->lpFirstNetwork,
		this->lpSecondNetwork);
	this->lpFirstNetworkCache = pCache->pNetworkCache(this->lpFirstNetwork);
}

}
