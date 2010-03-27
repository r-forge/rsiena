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
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructs a new network ministep.
 * @param[in] variableId the ID of the dependent variable to be changed
 * @param[in] ego the actor making the change
 * @param[in] alter the alter whose incoming tie is changed
 * @param[in] difference the amount of change
 * (-1,0,+1 for dichotomous variables)
 */
NetworkChange::NetworkChange(int variableId,
	int ego,
	int alter,
	int difference) : MiniStep(variableId, ego, difference)
{
	this->lalter = alter;
}


/**
 * Deallocates this ministep.
 */
NetworkChange::~NetworkChange()
{
}


/**
 * Changes the given network variable according to this ministep.
 */
void NetworkChange::makeChange(DependentVariable * pVariable)
{
	MiniStep::makeChange(pVariable);

	if (this->ego() != this->lalter &&
		this->difference() != 0)
	{
		NetworkVariable * pNetworkVariable =
			dynamic_cast<NetworkVariable *>(pVariable);
		int oldValue = pNetworkVariable->pNetwork()->tieValue(this->ego(),
			this->lalter);
		pNetworkVariable->pNetwork()->setTieValue(this->ego(),
			this->lalter,
			oldValue + this->difference());
	}
}


/**
 * Returns if this ministep is diagonal, namely, it does not change
 * the dependent variables.
 */
bool NetworkChange::diagonal() const
{
	return this->ego() == this->lalter || this->difference() == 0;
}

}
