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

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class GwespBFEffect : public NetworkEffect
{
public:
	GwespBFEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	double weight;
};

}

#endif /*GwespBFEffect_H_*/
