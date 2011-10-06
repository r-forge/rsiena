/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: DiffusionRateEffect.cpp
 *
 * Description: This file contains the implementation of the class
 * DiffusionRateEffect.
 *****************************************************************************/

#include "DiffusionRateEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/DiffusionEffectValueTable.h"
#include "network/IncidentTieIterator.h"
#include <R_ext/Print.h>

namespace siena
{

/**
 * Constructor.
 * @param[in] pVariable the network variable this effect depends on
 * @param[in] pBehaviorVariable the behavior variable this effect depends on
 * @param[in] type the type of this effect
 * @param[in] parameter the statistical parameter of this effect
 */

DiffusionRateEffect::DiffusionRateEffect(const NetworkVariable * pVariable,
	const BehaviorVariable * pBehaviorVariable,
	DiffusionRateEffectType type,
	double parameter)
{
	this->lpVariable = pVariable;
	this->lpBehaviorVariable = pBehaviorVariable;
	this->ltype = type;

	double possibleDegreeNumer = 1;
	double possibleDegreeDenom = 1;

	if (this->ltype == AVERAGE_EXPOSURE_RATE)
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = max(this->lpVariable->n(),
			this->lpVariable->m());
	}
	this->lpTable = new DiffusionEffectValueTable(possibleDegreeNumer,
		possibleDegreeDenom);
	this->lpTable->parameter(parameter);
}

/**
 * Destructor.
 */
DiffusionRateEffect::~DiffusionRateEffect()
{
	delete this->lpTable;
	this->lpTable = 0;
}


/**
 * Returns the contribution of this effect for the given actor.
 */
double DiffusionRateEffect::value(int i)
{
	Network * pNetwork = this->lpVariable->pNetwork();
	switch (this->ltype)
	{
		case AVERAGE_EXPOSURE_RATE:
		{
			int totalAlterValue = 0;
			if (pNetwork->outDegree(i) > 0)
			{
				for (IncidentTieIterator iter = pNetwork->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = this->lpBehaviorVariable->
						value(iter.actor());
					totalAlterValue += alterValue;
				}
			}

			if (totalAlterValue==0)
			{
				return 1;
			}
			if (totalAlterValue>0)
			{
				return this->lpTable->value(totalAlterValue,max(1,
						pNetwork->outDegree(i)));
			}
		}
	}

	throw new logic_error("Unexpected diffusion rate effect type");
}
/**
 * Stores the parameter for the diffusion rate effect.
 */
void DiffusionRateEffect::parameter(double parameterValue)
{
	this->lpTable->parameter(parameterValue);
}

/**
 * Returns the parameter for the diffusion rate effect.
 */
double DiffusionRateEffect::parameter() const
{
	return this->lpTable->parameter();
}


}
