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

#include <cmath>
#include "DiffusionRateEffect.h"
#include "utils/Utils.h"
#include "network/OneModeNetwork.h"
#include "model/variables/NetworkVariable.h"
#include "model/variables/BehaviorVariable.h"
#include "model/variables/DiffusionEffectValueTable.h"
#include "data/ConstantCovariate.h"
#include "data/ChangingCovariate.h"
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
	if (this->ltype == SUSCEPT_AVERAGE_INDEGREE_RATE)
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = max(this->lpVariable->n(),
			this->lpVariable->m());
	}
	if (this->ltype == TOTAL_EXPOSURE_RATE)
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = 1;
	}
	if (this->ltype == INFECTION_INDEGREE_RATE)
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = 1;
	}
	if (this->ltype == INFECTION_OUTDEGREE_RATE)
	{
	   	possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = 1;
	}
	this->lpTable = new DiffusionEffectValueTable(possibleDegreeNumer,
		possibleDegreeDenom);
	this->lpTable->parameter(parameter);
}

DiffusionRateEffect::DiffusionRateEffect(const NetworkVariable * pVariable,
	const BehaviorVariable * pBehaviorVariable,
	const ConstantCovariate * pConstantCovariate,
	const ChangingCovariate * pChangingCovariate,
	DiffusionRateEffectType type,
	double parameter)
{
	this->lpVariable = pVariable;
	this->lpBehaviorVariable = pBehaviorVariable;
	this->lpChangingCovariate = pChangingCovariate;
	this->lpConstantCovariate = pConstantCovariate;
	this->ltype = type;

	double possibleDegreeNumer = 1;
	double possibleDegreeDenom = 1;

	if (this->ltype == SUSCEPT_AVERAGE_COVARIATE_RATE)
	{
		possibleDegreeNumer = (this->lpBehaviorVariable->range()) *
			max(this->lpVariable->n(),
				this->lpVariable->m());
		possibleDegreeDenom = max(this->lpVariable->n(),
			this->lpVariable->m());
	}
	if (this->ltype == INFECTION_COVARIATE_RATE)
	{
		possibleDegreeNumer = 1;
		possibleDegreeDenom = 1;
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
double DiffusionRateEffect::value(int i, int period)
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
		case SUSCEPT_AVERAGE_INDEGREE_RATE:
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

			totalAlterValue *= pNetwork->inDegree(i);

			if (totalAlterValue==0)
			{
			    return 1;
			}
			else
			{
			    return this->lpTable->value(totalAlterValue,max(1,
						pNetwork->outDegree(i)));
			}
		}
		case TOTAL_EXPOSURE_RATE:
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
			else
			{
			    return this->lpTable->value(totalAlterValue, 1);
			}
		}
		case SUSCEPT_AVERAGE_COVARIATE_RATE:
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
			else
			{
				if (this->lpConstantCovariate)
				{
					return pow(this->lpTable->value(totalAlterValue,max(1,
								pNetwork->outDegree(i))),
						this->lpConstantCovariate->value(i));
				}
				else if (this->lpChangingCovariate)
				{
					return pow(this->lpTable->value(totalAlterValue,max(1,
								pNetwork->outDegree(i))),
						this->lpChangingCovariate->value(i,period));
				}
			   	else
					throw logic_error(
						"No individual covariate found.");
			}
		}
		case INFECTION_INDEGREE_RATE:
		{
			int totalAlterValue = 0;
			if (pNetwork->outDegree(i) > 0)
			{
				for (IncidentTieIterator iter = pNetwork->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = this->lpBehaviorVariable->
						value(iter.actor()) *
						pNetwork->inDegree(iter.actor());
					totalAlterValue += alterValue;
				}
			}

			if (totalAlterValue==0)
			{
			    return 1;
			}
			else
			{
			    return this->lpTable->value(totalAlterValue,1);
			}
		}
		case INFECTION_OUTDEGREE_RATE:
		{
			int totalAlterValue = 0;
			if (pNetwork->outDegree(i) > 0)
			{
				for (IncidentTieIterator iter = pNetwork->outTies(i);
					 iter.valid();
					 iter.next())
				{
					double alterValue = this->lpBehaviorVariable->
						value(iter.actor()) *
						pNetwork->outDegree(iter.actor());
					totalAlterValue += alterValue;
				}
			}

			if (totalAlterValue==0)
			{
			    return 1;
			}
			else
			{
				return this->lpTable->value(totalAlterValue,1);
			}
		}
		case INFECTION_COVARIATE_RATE:
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

					if (this->lpConstantCovariate)
					{
						alterValue *=
							this->lpConstantCovariate->value(iter.actor());
					}
					else if (this->lpChangingCovariate)
					{
						alterValue *=
							this->lpChangingCovariate->value(iter.actor(),
								period);
					}
					else
						throw logic_error(
							"No individual covariate found.");

					totalAlterValue += alterValue;
				}
			}
			if (totalAlterValue==0)
			{
			    return 1;
			}
			else
			{
			    return pow(this->lpTable->value(1,1), totalAlterValue);
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
