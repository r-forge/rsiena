/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.h
 *
 * Description: This file contains the definition of the
 * DiffusionRateEffect class.
 *****************************************************************************/

#ifndef DIFFUSIONRATEEFFECT_H_
#define DIFFUSIONRATEEFFECT_H_

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Enums
// ----------------------------------------------------------------------------

/**
 * Available types of diffusion rate effects.
 */
enum DiffusionRateEffectType
{
	AVERAGE_EXPOSURE_RATE	
};


// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class NetworkVariable;
class BehaviorVariable;
class DiffusionEffectValueTable;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

/**
 * Encapsulates the information necessary for calculating the contributions
 * of a diffusion rate effect. This includes the effect type,
 * the network variable and the behavior variable
 * the effect depends on, and the statistical parameter of the effect. 
 */
class DiffusionRateEffect
{
public:
	DiffusionRateEffect(const NetworkVariable * pVariable,
		const BehaviorVariable * pBehaviorVariable,
		DiffusionRateEffectType type,
		double parameter);
	virtual ~DiffusionRateEffect();

	double value(int i);
	void parameter(double parameterValue);
	double parameter() const;

private:
	// The network variable this effect depends on
	const NetworkVariable * lpVariable;

	// The behavior variable this effect depends on
	const BehaviorVariable * lpBehaviorVariable;

	// The type of the effect
	DiffusionRateEffectType ltype;

	// A table for efficient calculation of contributions. If two actors have
	// the same effect value, then the table ensures that we don't  
	// calculate the same contribution twice.

	DiffusionEffectValueTable * lpTable;

};

}

#endif /* DIFFUSIONRATEEFFECT_H_ */
