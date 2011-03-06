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

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class GwespFBEffect : public NetworkEffect
{
public:
	GwespFBEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	double weight;
};

}

#endif /*GwespFBEffect_H_*/
