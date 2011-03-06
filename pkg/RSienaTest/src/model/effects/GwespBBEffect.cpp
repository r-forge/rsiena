/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespBBEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespBBEffect.
 *****************************************************************************/

#include "GwespBBEffect.h"
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
GwespBBEffect::GwespBBEffect(
	const EffectInfo * pEffectInfo) : NetworkEffect(pEffectInfo)
{
	this->weight = -.01 * pEffectInfo->internalEffectParameter();
}


/**
 * Calculates the contribution of a tie flip to the given actor.
 */
double GwespBBEffect::calculateContribution(int alter) const
{
	return exp(this->weight * this->pReverseTwoPathTable()->get(alter) );
}


/**
 * The contribution of the tie from the implicit ego to the given alter
 * to the statistic. It is assumed that preprocessEgo(ego) has been
 * called before.
 */
double GwespBBEffect::tieStatistic(int alter)
{
	return exp(this->weight * this->pReverseTwoPathTable()->get(alter) ) / 6.0;
}

}
