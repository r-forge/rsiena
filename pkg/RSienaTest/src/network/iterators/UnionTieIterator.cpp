/*
 * UnionTieIterator.cpp
 *
 *  Created on: 30.07.2013
 *      Author: ortmann
 */

#include "UnionTieIterator.h"

namespace siena {

UnionTieIterator::UnionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		CombinedTieIterator(iter1, iter2) {
	init(1, 2);
}

UnionTieIterator::UnionTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2, int idIter1, int idIter2) :
		CombinedTieIterator(iter1, iter2) {
	init(idIter1, idIter2);
}

UnionTieIterator::~UnionTieIterator() {
}

void UnionTieIterator::next() {
	// at least one of the iterators has to be valid
	if (!valid()) {
		return;
	}
	// if both iterators are valid move to the next position
	if (lpIter1->valid() && lpIter2->valid()) {
		if (lpIter1->actor() < lpIter2->actor()) {
			lpIter1->next();
		} else if (lpIter1->actor() > lpIter2->actor()) {
			lpIter2->next();
		} else {
			lpIter1->next();
			lpIter2->next();
		}
	} else if (lpIter1->valid()) {
		lpIter1->next();
	} else {
		lpIter2->next();
	}
}

int UnionTieIterator::actor() const {
	// if both iterators are valid return the actor with the lower value
	if (lpIter1->valid() && lpIter2->valid()) {
		if (lpIter1->actor() <= lpIter2->actor()) {
			return lpIter1->actor();
		} else {
			return lpIter2->actor();
		}
		// else return the actor of the valid iterator
	} else if (lpIter1->valid()) {
		return lpIter1->actor();
	} else if (lpIter2->valid()) {
		return lpIter2->actor();
	}
	throw InvalidIteratorException();
}

bool UnionTieIterator::valid() const {
	return lpIter1->valid() || lpIter2->valid();
}

UnionTieIterator * UnionTieIterator::clone() const {
	return new UnionTieIterator(*this);
}

int UnionTieIterator::getInactiveIterID() {
	if (valid()) {
		if (!lpIter1->valid()) {
			return lIdIter1;
		}
		if (!lpIter2->valid() || lpIter1->actor() < lpIter2->actor()) {
			return lIdIter2;
		}
		return lIdIter1;
	}
	throw InvalidIteratorException();
}

int UnionTieIterator::getActiveIterID() {
	if (getInactiveIterID() == lIdIter1) {
		return lIdIter2;
	}
	return lIdIter1;
}

bool UnionTieIterator::isCommonNeighbor() const {
	return lpIter1->valid() && lpIter2->valid() && isCommon();
}

UnionTieIterator::UnionTieIterator(const UnionTieIterator& rhs) :
		CombinedTieIterator(rhs), //
		lIdIter1(rhs.lIdIter1), //
		lIdIter2(rhs.lIdIter2) {
}

void UnionTieIterator::init(int idIter1, int idIter2) {
	lIdIter1 = idIter1;
	lIdIter2 = idIter2;
}

} /* namespace siena */
