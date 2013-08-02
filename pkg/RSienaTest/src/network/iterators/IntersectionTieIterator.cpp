/*
 * IntersectionTieIterator.cpp
 *
 *  Created on: 30.07.2013
 *      Author: ortmann
 */

#include "IntersectionTieIterator.h"

namespace siena {

IntersectionTieIterator::IntersectionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		CombinedTieIterator(iter1, iter2) {
	if (valid() && !isCommon()) {
		skip();
	}
}

IntersectionTieIterator::~IntersectionTieIterator() {
}

int IntersectionTieIterator::actor() const {
	if (valid()) {
		return lpIter1->actor();
	}
	throw InvalidIteratorException();
}

void IntersectionTieIterator::next() {
	lpIter1->next();
	lpIter2->next();
	skip();
}

void IntersectionTieIterator::skip() {
	while (valid() && !isCommon()) {
		while (lpIter1->valid() && lpIter1->actor() < lpIter2->actor()) {
			lpIter1->next();
		}
		if (!lpIter1->valid()) {
			return;
		}
		while (lpIter2->valid() && lpIter2->actor() < lpIter1->actor()) {
			lpIter2->next();
		}
	}
}

bool IntersectionTieIterator::valid() const {
	return lpIter1->valid() && lpIter2->valid();
}

IntersectionTieIterator* IntersectionTieIterator::clone() const {
	return new IntersectionTieIterator(*this);
}

IntersectionTieIterator::IntersectionTieIterator(
		const IntersectionTieIterator& rhs) :
		CombinedTieIterator(rhs) {
	// note there is no need to skip as rhs is an IntersectionTieIterator itself
}

} /* namespace siena */
