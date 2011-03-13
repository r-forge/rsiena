/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFBEffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespFBEffect.
 *****************************************************************************/

#ifndef GwespFBEffect_H_
#define GwespFBEffect_H_

#include "model/effects/NetworkEffect.h"
#include "GwespEffect.h"

namespace siena
{

class GwespFBEffect : public GwespEffect
{
public:
	GwespFBEffect(const EffectInfo * pEffectInfo);
protected:
	virtual inline ConfigurationTable * pStatisticTable() const;
};

}

#endif /*GwespFBEffect_H_*/
