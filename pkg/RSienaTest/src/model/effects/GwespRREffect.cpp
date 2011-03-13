/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespRREffect.cpp
 *
 * Description: This file contains the implementation of the class
 * GwespRREffect.
 *****************************************************************************/

#include "GwespRREffect.h"

namespace siena
{

GwespRREffect::GwespRREffect(const EffectInfo * pEffectInfo) :
		GwespEffect(pEffectInfo)
{
	// Do nothing
}

/**
 * Uses the reciprocated-reciprocated table to calculate edgewise shared partners
 */
inline ConfigurationTable * GwespRREffect::pStatisticTable() const
{
	return this->pRRTable();
}



}
