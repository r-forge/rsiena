/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: LinearShapeEffect.cpp
 *
 * Description: This file contains the implementation of the
 * LinearShapeEffect class.
 *****************************************************************************/
#include <R.h>

#include "LinearShapeEffect.h"
#include "model/variables/BehaviorVariable.h"
#include "data/BehaviorLongitudinalData.h"

namespace siena
{

/**
 * Constructor.
 */
LinearShapeEffect::LinearShapeEffect(const EffectInfo * pEffectInfo) :
	BehaviorEffect(pEffectInfo)
{
}


/**
 * Calculates the change in the statistic corresponding to this effect if
 * the given actor would change his behavior by the given amount.
 */
double LinearShapeEffect::calculateChangeContribution(int actor,
	int difference) const
{
	return difference;
}

/**
 * Returns the statistic corresponding to this effect as part of
 * the evaluation function with respect to the given behavior variable.
 */
double LinearShapeEffect::evaluationStatistic(double * currentValues) const
{
	double statistic = 0;
	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		statistic += currentValues[i];
	}

	return statistic;
}


/**
 * Returns the statistic corresponding to this effect as part of
 * the endowment function with respect to the initial values of a
 * behavior variable and the current values.
 */
double LinearShapeEffect::endowmentStatistic(const int * difference,
	double * currentValues) const
{
	double statistic = 0;

	int n = this->n();

	for (int i = 0; i < n; i++)
	{
		if (difference[i] > 0)
		{
			statistic += currentValues[i]  ;
		}
	}

	return statistic;
}

}