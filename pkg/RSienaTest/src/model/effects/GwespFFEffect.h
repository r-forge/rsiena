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

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class GwespFFEffect : public NetworkEffect
{
public:
	GwespFFEffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	double weight;
};

}

#endif /*GwespFFEffect_H_*/
