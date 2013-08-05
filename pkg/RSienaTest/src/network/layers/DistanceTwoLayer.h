/*
 * DistanceTwoLayer.h
 *
 *  Created on: 05.08.2013
 *      Author: ortmann
 */

#ifndef DISTANCETWOLAYER_H_
#define DISTANCETWOLAYER_H_

#include <map>

#include "NetworkLayer.h"
#include "../Network.h"

namespace siena {

class DistanceTwoLayer: public NetworkLayer {
public:
	DistanceTwoLayer(const Network& rNetwork);
	virtual ~DistanceTwoLayer();

	void initializeOneMode(const Network& rNetwork);

	void initializeTwoMode(const Network& rNetwork);

	void onTieIntroductionEvent(const Network& rNetwork, const int ego,
			const int alter);
	void onTieWithdrawalEvent(const Network& rNetwork, const int ego,
			const int alter);

	void onNetworkClearEvent(const Network& rNetwork);

	IncidentTieIterator getDistanceTwoNeighbors(int ego) const;
private:
	void updateSingleTieValue(int ego, int alter, int val);

	void modifyTieValue(int ego, int alter, int val);

	void modify2PathCountOneMode(const Network& rNetwork, int ego, int alter,
			int val);

	void modify2PathCountTwoMode(const Network& rNetwork, int ego, int alter,
			int val);

	std::map<int, int>* lpAdjacencies;
};

} /* namespace siena */
#endif /* DISTANCETWOLAYER_H_ */
