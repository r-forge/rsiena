/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFBEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespFBEffect.
 *****************************************************************************/

#include "GwespFBEffect.h"

namespace siena
{

GwespFBEffect::GwespFBEffect(const EffectInfo * pEffectInfo) :
		GwespEffect(pEffectInfo)
{
	// Do nothing
}
/**
 * Uses the in star table to calculate edgewise shared partners
 */
inline ConfigurationTable * GwespFBEffect::pStatisticTable() const
{
	return this->pInStarTable();
}



}
