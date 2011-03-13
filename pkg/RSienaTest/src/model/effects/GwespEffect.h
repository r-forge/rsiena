/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespEffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespEffect.
 *****************************************************************************/

#ifndef GwespEffect_H_
#define GwespEffect_H_

#define MAX_STATISTIC 1000

#include "model/effects/NetworkEffect.h"
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "model/variables/NetworkVariable.h"
#include "model/tables/ConfigurationTable.h"
#include "model/EffectInfo.h"
#include <math.h>

namespace siena
{

/**
 * This abstract class defines the Gwesp effect.
 * All Gwesp effects implement
 * this class by e.g. pointing to the correct ConfigurationTable
 * on the inline function pStatisticTable().
 */
class GwespEffect : public NetworkEffect
{
public:
	GwespEffect(const EffectInfo * pEffectInfo);
	double calculateContribution(int alter) const;

protected:
	virtual inline ConfigurationTable * pStatisticTable() const = 0;
	double tieStatistic(int alter);
	double weight;
	double cumulativeWeight[MAX_STATISTIC];
};

}

#endif /*GwespEffect_H_*/
