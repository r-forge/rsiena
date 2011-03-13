/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespFFEffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespFFEffect.
 *****************************************************************************/

#ifndef GwespFFEffect_H_
#define GwespFFEffect_H_

#include "model/effects/NetworkEffect.h"
#include "GwespEffect.h"

namespace siena
{

class GwespFFEffect : public GwespEffect
{
public:
	GwespFFEffect(const EffectInfo * pEffectInfo);
protected:
	virtual inline ConfigurationTable * pStatisticTable() const;
};

}

#endif /*GwespFFEffect_H_*/
