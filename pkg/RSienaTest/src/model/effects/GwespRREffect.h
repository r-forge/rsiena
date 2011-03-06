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

namespace siena
{

/**
 * This class defines the transitive triads effect. It is a version of the
 * transitive triplets effect for symmetric networks.
 */
class GwespRREffect : public NetworkEffect
{
public:
	GwespRREffect(const EffectInfo * pEffectInfo);

	virtual double calculateContribution(int alter) const;

protected:
	virtual double tieStatistic(int alter);
	double weight;
};

}

#endif /*GwespRREffect_H_*/
