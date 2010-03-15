/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: Chain.cpp
 *
 * Description: This file contains the implementation of the class Chain.
 *****************************************************************************/

#include <vector>
#include "Chain.h"
#include "utils/Random.h"
#include "model/ml/MiniStep.h"
#include "model/ml/BehaviorChange.h"
#include "model/ml/NetworkChange.h"
#include "model/ml/MLSimulation.h"
#include "network/Network.h"
#include "network/IncidentTieIterator.h"
#include "data/Data.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/NetworkLongitudinalData.h"

namespace siena
{

/**
 * Creates an empty chain.
 */
Chain::Chain(Data * pData)
{
	this->lpFirst = new MiniStep(0, 0, 0);
	this->lpLast = new MiniStep(0, 0, 0);

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->lpData = pData;
	this->lperiod = -1;

	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->lmu = 0;
	this->lsigma2 = 0;
}


/**
 * Deallocates this chain.
 */
Chain::~Chain()
{
	this->clear();

	delete this->lpFirst;
	delete this->lpLast;
	this->lpFirst = 0;
	this->lpLast = 0;
	this->lpData = 0;

	this->lminiSteps.clear();
	this->ldiagonalMiniSteps.clear();
}


/**
 * Removes all ministeps from this chain except for the dummy ministeps
 * at the ends of the chain.
 */
void Chain::clear()
{
	MiniStep * pMiniStep = this->lpFirst->pNext();

	while (pMiniStep != this->lpLast)
	{
		MiniStep * pNextMiniStep = pMiniStep->pNext();
		delete pMiniStep;
		pMiniStep = pNextMiniStep;
	}

	this->lpFirst->pNext(this->lpLast);
	this->lpLast->pPrevious(this->lpFirst);

	this->lminiSteps.clear();
	this->lminiSteps.push_back(this->lpLast);
	this->lpLast->index(0);

	this->ldiagonalMiniSteps.clear();

	this->lmu = 0;
	this->lsigma2 = 0;
}


/**
 * Inserts a new ministep before a ministep that already belongs to this chain.
 */
void Chain::insertBefore(MiniStep * pNewMiniStep, MiniStep * pExistingMiniStep)
{
	MiniStep * pPreviousMiniStep = pExistingMiniStep->pPrevious();

	pNewMiniStep->pChain(this);

	pNewMiniStep->pPrevious(pPreviousMiniStep);
	pPreviousMiniStep->pNext(pNewMiniStep);

	pNewMiniStep->pNext(pExistingMiniStep);
	pExistingMiniStep->pPrevious(pNewMiniStep);

	if (pNewMiniStep->diagonal())
	{
		this->ldiagonalMiniSteps.push_back(pNewMiniStep);
		pNewMiniStep->diagonalIndex(this->ldiagonalMiniSteps.size() - 1);
	}

	this->lminiSteps.push_back(pNewMiniStep);
	pNewMiniStep->index(this->lminiSteps.size() - 1);

	double rr = pNewMiniStep->reciprocalRate();
	this->lmu += rr;
	this->lsigma2 += rr * rr;
}


/**
 * Removes the given ministep from this chain without deleting the ministep.
 */
void Chain::remove(MiniStep * pMiniStep)
{
	pMiniStep->pPrevious()->pNext(pMiniStep->pNext());
	pMiniStep->pNext()->pPrevious(pMiniStep->pPrevious());
	pMiniStep->pNext(0);
	pMiniStep->pPrevious(0);

	if (pMiniStep->diagonal())
	{
		MiniStep * pLastMiniStep =
			this->ldiagonalMiniSteps[this->ldiagonalMiniSteps.size() - 1];
		this->ldiagonalMiniSteps[pMiniStep->diagonalIndex()] = pLastMiniStep;
		pLastMiniStep->diagonalIndex(pMiniStep->diagonalIndex());
		this->ldiagonalMiniSteps.pop_back();
		pMiniStep->diagonalIndex(-1);
	}

	MiniStep * pLastMiniStep = this->lminiSteps[this->lminiSteps.size() - 1];
	this->lminiSteps[pMiniStep->index()] = pLastMiniStep;
	pLastMiniStep->index(pMiniStep->index());
	this->lminiSteps.pop_back();
	pMiniStep->index(-1);

	double rr = pMiniStep->reciprocalRate();
	this->lmu -= rr;
	this->lsigma2 -= rr * rr;

	pMiniStep->pChain(0);
}


/**
 * Generates a random chain connecting the start and end observations of the
 * given data object for the given period. The chain is simple in the sense
 * that no two ministeps cancel each other out.
 */
void Chain::connect(int period)
{
	this->clear();
	this->lperiod = period;
	vector<MiniStep *> miniSteps;

	// Create the required ministeps

	for (unsigned variableIndex = 0;
		variableIndex < this->lpData->rDependentVariableData().size();
		variableIndex++)
	{
		LongitudinalData * pVariableData =
			this->lpData->rDependentVariableData()[variableIndex];
		NetworkLongitudinalData * pNetworkData =
			dynamic_cast<NetworkLongitudinalData *>(pVariableData);
		BehaviorLongitudinalData * pBehaviorData =
			dynamic_cast<BehaviorLongitudinalData *>(pVariableData);

		if (pNetworkData)
		{
			const Network * pNetwork1 = pNetworkData->pNetwork(period);
			const Network * pNetwork2 = pNetworkData->pNetwork(period + 1);

			for (int i = 0; i < pNetwork1->n(); i++)
			{
				IncidentTieIterator iter1 = pNetwork1->outTies(i);
				IncidentTieIterator iter2 = pNetwork2->outTies(i);

				while (iter1.valid() || iter2.valid())
				{
					if (iter1.valid() &&
						(!iter2.valid() || iter1.actor() < iter2.actor()))
					{
						miniSteps.push_back(
							new NetworkChange(i,
								iter1.actor(),
								pNetworkData->name(),
								-1));
						iter1.next();
					}
					else if (iter2.valid() &&
						(!iter1.valid() || iter2.actor() < iter1.actor()))
					{
						miniSteps.push_back(
							new NetworkChange(i,
								iter2.actor(),
								pNetworkData->name(),
								1));
						iter2.next();
					}
					else
					{
						iter1.next();
						iter2.next();
					}
				}
			}
		}
		else if (pBehaviorData)
		{
			for (int i = 0; i < pBehaviorData->n(); i++)
			{
				int delta = pBehaviorData->value(period + 1, i) -
					pBehaviorData->value(period, i);
				int singleChange = 1;

				if (delta < 0)
				{
					delta = -delta;
					singleChange = -1;
				}

				for (int j = 0; j < delta; j++)
				{
					miniSteps.push_back(
						new BehaviorChange(i,
							pBehaviorData->name(),
							singleChange));
				}
			}
		}
	}

	// Randomize the ministeps

	for (unsigned i = 1; i < miniSteps.size(); i++)
	{
		int j = nextInt(i + 1);
		MiniStep * pTempMiniStep = miniSteps[i];
		miniSteps[i] = miniSteps[j];
		miniSteps[j] = pTempMiniStep;
	}

	// And finally add the ministeps to this chain

	for (unsigned i = 0; i < miniSteps.size(); i++)
	{
		this->insertBefore(miniSteps[i], this->lpLast);
	}
}


/**
 * Performs the necessary updates when the reciprocal rate of the
 * given ministep changes to the given value.
 */
void Chain::onReciprocalRateChange(const MiniStep * pMiniStep, double newValue)
{
	double oldValue = pMiniStep->reciprocalRate();

	this->lmu -= oldValue;
	this->lsigma2 -= oldValue * oldValue;

	this->lmu += newValue;
	this->lsigma2 += newValue * newValue;
}


/**
 * Returns the period whose end observations are connected by this chain.
 */
int Chain::period() const
{
	return this->lperiod;
}


/**
 * Returns the first (dummy) ministep in this chain.
 */
MiniStep * Chain::pFirst() const
{
	return this->lpFirst;
}


/**
 * Returns the last (dummy) ministep in this chain.
 */
MiniStep * Chain::pLast() const
{
	return this->lpLast;
}


/**
 * Returns the number of ministeps of this chain excluding the first,
 * but including the last dummy ministep.
 */
int Chain::ministepCount() const
{
	return this->lminiSteps.size();
}


/**
 * Returns the number of diagonal ministeps in this chain.
 */
int Chain::diagonalMinistepCount() const
{
	return this->ldiagonalMiniSteps.size();
}


/**
 * Returns the sum of reciprocal rates over all non-dummy ministeps.
 */
double Chain::mu() const
{
	return this->lmu;
}


/**
 * Returns the sum of squared reciprocal rates over all non-dummy ministeps.
 */
double Chain::sigma2() const
{
	return this->lsigma2;
}


// ----------------------------------------------------------------------------
// Section: Random draws
// ----------------------------------------------------------------------------

/**
 * Returns a random ministep excluding the first dummy ministep.
 */
MiniStep * Chain::randomMiniStep() const
{
	return this->lminiSteps[nextInt(this->lminiSteps.size())];
}


/**
 * Returns a random diagonal ministep.
 */
MiniStep * Chain::randomDiagonalMiniStep() const
{
	return this->ldiagonalMiniSteps[nextInt(this->ldiagonalMiniSteps.size())];
}

}
