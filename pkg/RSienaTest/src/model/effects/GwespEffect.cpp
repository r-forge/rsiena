/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespEffect.
 *****************************************************************************/

#include "GwespEffect.h"

namespace siena
{

/**
 * Constructor.
 */
GwespEffect::GwespEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->weight = -.01 * pEffectInfo->internalEffectParameter();
	cumulativeWeight[0] = 0;
	for (int i=1; i<MAX_STATISTIC; i++) {
		cumulativeWeight[i] = cumulativeWeight[i-1] + exp(this->weight * (double)i);
	}
	// And normalize
	for (int i=1; i<MAX_STATISTIC; i++) {
		cumulativeWeight[i] /= cumulativeWeight[MAX_STATISTIC-1];
	}
}

/**
 * Calculates the change statistic for any Gwesp effect
 * @param[in] alter The alter under consideration
 */
double GwespEffect::calculateContribution(int alter) const
{
	int statistic = this->pStatisticTable()->get(alter);
	if (statistic >= MAX_STATISTIC) {
		return 1.0;
	} else {
		return cumulativeWeight[statistic];
	}
}

/**
 * Calculates the tie statistic for any Gwesp effect
 * @param[in] alter The which alter to calculate edgewise shared partners
 */
double GwespEffect::tieStatistic(int alter)
{
	int statistic = this->pStatisticTable()->get(alter);
	if (statistic >= MAX_STATISTIC) {
		return 1.0;
	} else {
		return cumulativeWeight[statistic];
	}
}

}
