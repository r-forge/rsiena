/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: CovariateDiffEffect.cpp
 *
 * Description: This file contains the implementation of the
 * CovariateDiffEffect class.
 *****************************************************************************/

#include "CovariateDiffEffect.h"
#include "network/Network.h"
#include "model/variables/NetworkVariable.h"

namespace siena
{

/**
 * Constructor.
 * @param[in] pEffectInfo the effect descriptor
 * @param[in] squared indicates if the covariate values must be squared
 */
CovariateDiffEffect::CovariateDiffEffect(const EffectInfo * pEffectInfo,
	bool squared) :
		CovariateDependentNetworkEffect(pEffectInfo)
{
	this->lsquared = squared;
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double CovariateDiffEffect::calculateContribution(int alter) const
{
	double change = this->value(alter) - this->value(this->ego());

	if (this->lsquared)
	{
		change *= change;
	}

	return change;
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double CovariateDiffEffect::tieStatistic(int alter)
{
	double statistic = 0;

	if (!(this->missing(alter) || this->missing(this->ego())))
	{
		statistic = this->value(alter) - this->value(this->ego());

		if (this->lsquared)
		{
			statistic *= statistic;
		}
	}

	return statistic;
}

}
