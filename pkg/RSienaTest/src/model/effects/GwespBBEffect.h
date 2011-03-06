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

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class GwespBBEffect : public NetworkEffect
{
public:
	GwespBBEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	double weight;
};

}

#endif /*GwespBBEffect_H_*/
