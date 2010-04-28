/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkChange.cpp
 *
 * Description: This file contains the implementation of the class
 * NetworkChange.
 *****************************************************************************/

#include "NetworkChange.h"
#include "network/Network.h"
#include "data/NetworkLongitudinalData.h"
#include "model/variables/NetworkVariable.h"
#include "model/ml/Option.h"

namespace siena
{

/**
 * Constructs a new network ministep.
 * @param[in] pData the longitudinal data object for the
 * corresponding network variable
 * @param[in] ego the actor making the change
 * @param[in] alter the alter whose incoming tie is changed
 */
NetworkChange::NetworkChange(NetworkLongitudinalData * pData,
	int ego,
	int alter) : MiniStep(pData, ego)
{
	this->lpData = pData;
	this->lalter = alter;
	this->pOption(new Option(pData->id(), ego, alter));
}


/**
 * Deallocates this ministep.
 */
NetworkChange::~NetworkChange()
{
}


/**
 * Returns if this ministep is changing a network variable.
 */
bool NetworkChange::networkMiniStep() const
{
	return true;
}


/**
 * Changes the given network variable according to this ministep.
 */
void NetworkChange::makeChange(DependentVariable * pVariable)
{
	MiniStep::makeChange(pVariable);

	if (this->ego() != this->lalter)
	{
		NetworkVariable * pNetworkVariable =
			dynamic_cast<NetworkVariable *>(pVariable);
		int oldValue = pNetworkVariable->pNetwork()->tieValue(this->ego(),
			this->lalter);
		pNetworkVariable->pNetwork()->setTieValue(this->ego(),
			this->lalter,
			1 - oldValue);
	}
}


/**
 * Returns if this ministep is diagonal, namely, it does not change
 * the dependent variables.
 */
bool NetworkChange::diagonal() const
{
	return this->ego() == this->lalter;
}


/**
 * Returns if the observed data for this ministep is missing at
 * either end of the given period.
 */
bool NetworkChange::missing(int period) const
{
	return this->lpData->missing(this->ego(), this->lalter, period) ||
		this->lpData->missing(this->ego(), this->lalter, period + 1);
}


/**
 * Returns a new ministep that reverses the effect of this ministep.
 */
MiniStep * NetworkChange::createReverseMiniStep() const
{
	return new NetworkChange(this->lpData,
		this->ego(),
		this->lalter);
}

}
