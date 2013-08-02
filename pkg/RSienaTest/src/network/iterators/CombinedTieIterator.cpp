/*
 * CombinedTieIterator.cpp
 *
 *  Created on: 30.07.2013
 *      Author: ortmann
 */

#include "CombinedTieIterator.h"

namespace siena {

CombinedTieIterator::CombinedTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		ITieIterator(), //
		lpIter1(iter1.clone()), //
		lpIter2(iter2.clone()) {
}

CombinedTieIterator::~CombinedTieIterator() {
	delete lpIter1;
	delete lpIter2;
}

CombinedTieIterator::CombinedTieIterator(const CombinedTieIterator& rhs) :
		ITieIterator(rhs), //
		lpIter1(rhs.lpIter1->clone()), //
		lpIter2(rhs.lpIter2->clone()) {
}

} /* namespace siena */
