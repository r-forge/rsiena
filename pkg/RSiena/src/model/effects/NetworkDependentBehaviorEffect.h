/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkDependentBehaviorEffect.h
 *
 * Description: This file contains the definition of the
 * NetworkDependentBehaviorEffect class.
 *****************************************************************************/

#ifndef NETWORKDEPENDENTBEHAVIOREFFECT_H_
#define NETWORKDEPENDENTBEHAVIOREFFECT_H_

#include "BehaviorEffect.h"
#include "model/variables/BehaviorVariable.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Network;


// ----------------------------------------------------------------------------
// Section: NetworkDependentBehaviorEffect class
// ----------------------------------------------------------------------------

/**
 * The base class for all behavior effects depending on some network variable.
 */
class NetworkDependentBehaviorEffect : public BehaviorEffect
{
public:
	NetworkDependentBehaviorEffect(const EffectInfo * pEffectInfo);
	virtual ~NetworkDependentBehaviorEffect();

	virtual void initialize(const Data * pData,
		State * pState,
		int period,
		Cache * pCache);

protected:
	inline const Network * pNetwork() const;
	double totalAlterValue(int i) const;
	double totalInAlterValue(int i) const;
	virtual void preprocessEgo(int ego);

private:
	// The network this effect is interacting with
	const Network * lpNetwork;
	// total out- and in-alter values
	double * ltotalAlterValues;
	double * ltotalInAlterValues;
};


// ----------------------------------------------------------------------------
// Section: Inline methods
// ----------------------------------------------------------------------------

/**
 * Returns the network this effect is interacting with.
 */
const Network * NetworkDependentBehaviorEffect::pNetwork()
	const
{
	return this->lpNetwork;
}

}

#endif /*NETWORKDEPENDENTBEHAVIOREFFECT_H_*/
