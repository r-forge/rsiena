/*
 * DistanceTwoLayer.cpp
 *
 *  Created on: 05.08.2013
 *      Author: ortmann
 */

#include "DistanceTwoLayer.h"

#include <vector>
#include <math.h>

#include "../IncidentTieIterator.h"
#include "../iterators/UnionTieIterator.h"

namespace siena {

using namespace std;

typedef map<int, int> TieMap;

DistanceTwoLayer::DistanceTwoLayer(const Network& rNetwork) :
		NetworkLayer(), //
		lpAdjacencies(new map<int, int> [rNetwork.n()]) {
	if (rNetwork.isOneMode()) {
		initializeOneMode(rNetwork);
	} else {
		initializeTwoMode(rNetwork);
	}
}

DistanceTwoLayer::~DistanceTwoLayer() {
	delete[] lpAdjacencies;
}

void DistanceTwoLayer::initializeOneMode(const Network& rNetwork) {
	for (int i = 0; i < rNetwork.n(); ++i) {
		std::vector<int> neighAtDistOne;
		// avoid the time to copy
		neighAtDistOne.reserve(rNetwork.outDegree(i) + rNetwork.inDegree(i));
		// we could do this all with UnionIterators but it is much slower
		for (UnionTieIterator iter(rNetwork.inTies(i), rNetwork.outTies(i));
				iter.valid(); iter.next()) {
			// take care of loops
			if (iter.actor() != i) {
				neighAtDistOne.push_back(iter.actor());
			}
		}
		vector<int>::const_iterator iterEnd = neighAtDistOne.end();
		for (vector<int>::const_iterator outerIter = neighAtDistOne.begin();
				outerIter != iterEnd; ++outerIter) {
			int ego = *outerIter;
			for (vector<int>::const_iterator innerIter = outerIter + 1;
					innerIter != iterEnd; ++innerIter) {
				modifyTieValue(ego, *innerIter, 1);
			}
		}
	}
}

void DistanceTwoLayer::initializeTwoMode(const Network& rNetwork) {
	// this is a two mode network so we do not need to check for loops
	for (int i = 0; i < rNetwork.m(); ++i) {
		for (IncidentTieIterator outerIter = rNetwork.inTies(i);
				outerIter.valid(); outerIter.next()) {
			int outerActor = outerIter.actor();
			// copy the iterator
			IncidentTieIterator innerIter(outerIter);
			// move to the next position
			innerIter.next();
			for (; innerIter.valid(); innerIter.next()) {
				modifyTieValue(outerActor, innerIter.actor(), 1);
			}
		}
	}
}

void DistanceTwoLayer::modifyTieValue(int ego, int alter, int val) {
	updateSingleTieValue(ego, alter, val);
	updateSingleTieValue(alter, ego, val);
}

void DistanceTwoLayer::modify2PathCountOneMode(const Network& rNetwork, int ego,
		int alter, int val) {
	// if it is a loop or the edge (alter,ego) exists we have nothing to do
	if (ego == alter || rNetwork.hasEdge(alter, ego)) {
		return;
	}
	for (UnionTieIterator iter(
			UnionTieIterator(rNetwork.outTies(ego), rNetwork.inTies(ego)),
			UnionTieIterator(rNetwork.outTies(alter), rNetwork.inTies(alter)),
			ego, alter); iter.valid(); iter.next()) {
		int curNeighbor = iter.actor();
		// check whether the current neighbor is ego or alter itself
		if (curNeighbor != ego && curNeighbor != alter) {
			// if it is a common neighbor it creates a new 2-path with both
			if (iter.isCommonNeighbor()) {
				modifyTieValue(curNeighbor, ego, val);
				modifyTieValue(curNeighbor, alter, val);
			} else {
				// get the inactive iter id (ego or alter)
				int inactiveIterID = iter.getInactiveIterID();
				modifyTieValue(curNeighbor, inactiveIterID, val);
			}
		}
	}
}

void DistanceTwoLayer::modify2PathCountTwoMode(const Network& rNetwork, int ego,
		int alter, int val) {
	for (IncidentTieIterator iter = rNetwork.inTies(alter); iter.valid();
			iter.next()) {
		if (iter.actor() != ego) {
			modifyTieValue(ego, iter.actor(), val);
		}
	}
}

void DistanceTwoLayer::onTieIntroductionEvent(const Network& rNetwork,
		const int ego, const int alter) {
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, 1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

void DistanceTwoLayer::onTieWithdrawalEvent(const Network& rNetwork,
		const int ego, const int alter) {
	if (rNetwork.isOneMode()) {
		modify2PathCountOneMode(rNetwork, ego, alter, -1);
	} else {
		modify2PathCountTwoMode(rNetwork, ego, alter, -1);
	}
}

void DistanceTwoLayer::onNetworkClearEvent(const Network& rNetwork) {
	for (int i = 0; i < rNetwork.n(); ++i) {
		lpAdjacencies[i].clear();
	}
}

IncidentTieIterator DistanceTwoLayer::getDistanceTwoNeighbors(int ego) const {
	return IncidentTieIterator(lpAdjacencies[ego]);
}

void DistanceTwoLayer::updateSingleTieValue(int ego, int alter, int val) {
	TieMap& egoMap = lpAdjacencies[ego];
	TieMap::iterator iter = egoMap.lower_bound(alter);
	// we found the element
	if (iter != egoMap.end() && !egoMap.key_comp()(alter, iter->first)) {
		int newVal = iter->second + val;
		if (newVal) {
			iter->second = newVal;
		} else {
			egoMap.erase(iter);
		}
	} else {
		egoMap.insert(iter, TieMap::value_type(alter, val));
	}
}

} /* namespace siena */
