/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespBFEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespBFEffect.
 *****************************************************************************/

#include "GwespBFEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"
#include "model/EffectInfo.h"
#include <math.h>

namespace siena
{

/**
 * Constructor.
 */
GwespBFEffect::GwespBFEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->weight = -.01 * pEffectInfo->internalEffectParameter();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double GwespBFEffect::calculateContribution(int alter) const
{
	return exp(this->weight * this->pOutStarTable()->get(alter) );
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double GwespBFEffect::tieStatistic(int alter)
{
	return exp(this->weight * this->pOutStarTable()->get(alter) ) / 6.0;
}

}
