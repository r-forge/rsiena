/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: GwespRREffect.h
 *
 * Description: This file contains the declaration of the class
 * GwespRREffect.
 *****************************************************************************/

#ifndef GwespRREffect_H_
#define GwespRREffect_H_

#include "model/effects/NetworkEffect.h"
#include "GwespEffect.h"

namespace siena
{

class GwespRREffect : public GwespEffect
{
public:
	GwespRREffect(const EffectInfo * pEffectInfo);
protected:
	virtual inline ConfigurationTable * pStatisticTable() const;
};

}

#endif /*GwespRREffect_H_*/
