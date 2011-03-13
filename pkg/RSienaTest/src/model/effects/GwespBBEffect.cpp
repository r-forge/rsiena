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

namespace siena
{

GwespBBEffect::GwespBBEffect(const EffectInfo * pEffectInfo) :
		GwespEffect(pEffectInfo)
{
	// Do nothing
}
/**
 * Uses the reverse two path table to calculate edgewise shared partners
 */
inline ConfigurationTable * GwespBBEffect::pStatisticTable() const
{
	return this->pReverseTwoPathTable();
}



}
