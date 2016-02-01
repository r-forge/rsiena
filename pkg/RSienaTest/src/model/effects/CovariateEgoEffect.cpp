/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateEgoEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateEgoEffect class.
 *****************************************************************************/

#include "CovariateEgoEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 *
 * @param[in] pEffectInfo the effect descriptor
 */
CovariateEgoEffect::CovariateEgoEffect(const EffectInfo * pEffectInfo) :
	CovariateDependentNetworkEffect(pEffectInfo) {
}

/**
 * Constructor.
 *
 * @param[in] pEffectInfo the effect descriptor
 * @param simulatedState If `true` the value() function uses the simulated
 *        state, if any or the value at the end of the period.
 */
CovariateEgoEffect::CovariateEgoEffect(const EffectInfo * pEffectInfo,
		const bool simulatedState) :
	CovariateDependentNetworkEffect(pEffectInfo, simulatedState) {
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateEgoEffect::calculateContribution(int alter) const
{
	return this->value(this->ego());
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateEgoEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!this->missing(this->ego()))
	{
		statistic = this->value(this->ego());
	}

	return statistic;
}


/**
 * Returns if this effect is an ego effect.
 */
bool CovariateEgoEffect::egoEffect() const
{
	return true;
}

}
