/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFFEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespFFEffect.
 *****************************************************************************/

#include "GwespFFEffect.h"

namespace siena
{

GwespFFEffect::GwespFFEffect(const EffectInfo * pEffectInfo) :
		GwespEffect(pEffectInfo)
{
	// Do nothing
}
/**
 * Uses the two path table to calculate edgewise shared partners
 */
inline ConfigurationTable * GwespFFEffect::pStatisticTable() const
{
	return this->pTwoPathTable();
}



}
