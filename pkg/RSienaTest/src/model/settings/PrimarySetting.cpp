/*
 * PrimarySetting.cpp
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#include "PrimarySetting.h"

#include "../../network/Network.h"
#include "../../network/OneModeNetwork.h"
#include "../../network/iterators/ITieIterator.h"
#include "../../network/iterators/UnionTieIterator.h"
#include "../../network/iterators/FilteredIterator.h"
#include "../../network/iterators/SingleIterator.h"
#include "../../network/IncidentTieIterator.h"
#include "../../logger/Logger.h"

using namespace siena::logger;

namespace siena {

PrimarySetting::PrimarySetting() :
		GeneralSetting(), //
		rDistTwoLayer(), //
		lpNetwork(0), //
		lpiter(0) {
}

PrimarySetting::~PrimarySetting() {
	if (lpiter != 0) {
		delete lpiter;
	}
}

void PrimarySetting::initSetting(Network* const lpNetwork) {
	lpNetwork->addNetworkChangeListener(&rDistTwoLayer);
	this->lpNetwork = lpNetwork;
}

void PrimarySetting::terminateSetting(Network* const lpNetwork) {
	lpNetwork->removeNetworkChangeListener(&rDistTwoLayer);
	rDistTwoLayer.clear(lpNetwork->n());
	this->lpNetwork = 0;
}

void PrimarySetting::initSetting() {
	if (lpiter == 0) {
		IncidentTieIterator iter1 = lpNetwork->inTies(ego());
		IncidentTieIterator iter2 = lpNetwork->outTies(ego());
		UnionTieIterator uIter1(iter1, iter2);
		iter1 = rDistTwoLayer.getDistanceTwoNeighbors(ego());
		SingleIterator egoIter(ego());
		UnionTieIterator uIter2(iter1, egoIter);
		lpiter = new UnionTieIterator(uIter1, uIter2);
	} else {
		LOGS(Priority::ERROR)<< "setting has not been terminated\n";
		throw "setting has not been terminated";
	}
}

void PrimarySetting::terminateSetting() {
	if (lpiter != 0) {
	delete lpiter;
	lpiter=0;
	GeneralSetting::terminateSetting();
	} else {
		LOGS(Priority::ERROR)<< "setting has not been initialized\n";
		throw "setting has not been initialized";
	}
}

ITieIterator* PrimarySetting::getSteps() {
	if (lpiter != 0) {
		return lpiter->clone();
	}
	LOGS(Priority::ERROR)<< "setting has not been initialized\n";
	throw "setting has not been initialized";
}

int PrimarySetting::getSize() {
	return lpiter->size();
}

} /* namespace siena */
