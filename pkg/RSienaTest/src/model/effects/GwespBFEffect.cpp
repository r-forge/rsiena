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

namespace siena
{

GwespBFEffect::GwespBFEffect(const EffectInfo * pEffectInfo) :
		GwespEffect(pEffectInfo)
{
	// Do nothing
}
/**
 * Uses the out star table to calculate edgewise shared partners
 */
inline ConfigurationTable * GwespBFEffect::pStatisticTable() const
{
	return this->pOutStarTable();
}



}
