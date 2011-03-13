/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespBBEffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespBBEffect.
 *****************************************************************************/

#ifndef GwespBBEffect_H_
#define GwespBBEffect_H_

#include "model/effects/NetworkEffect.h"
#include "GwespEffect.h"

namespace siena
{

class GwespBBEffect : public GwespEffect
{
public:
	GwespBBEffect(const EffectInfo * pEffectInfo);
protected:
	virtual inline ConfigurationTable * pStatisticTable() const;
};

}

#endif /*GwespBBEffect_H_*/
