/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespBFEffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespBFEffect.
 *****************************************************************************/

#ifndef GwespBFEffect_H_
#define GwespBFEffect_H_

#include "model/effects/NetworkEffect.h"
#include "GwespEffect.h"

namespace siena
{

class GwespBFEffect : public GwespEffect
{
public:
	GwespBFEffect(const EffectInfo * pEffectInfo);
protected:
	virtual inline ConfigurationTable * pStatisticTable() const;
};

}

#endif /*GwespBFEffect_H_*/
