/*
 * SymDiffTieIterator.cpp
 *
 *  Created on: 30.07.2013
 *      Author: ortmann
 */

#include "SymDiffTieIterator.h"

namespace siena {

SymDiffTieIterator::SymDiffTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2) :
		UnionTieIterator(iter1, iter2) {
	init();
}

SymDiffTieIterator::SymDiffTieIterator(const ITieIterator& iter1,
		const ITieIterator& iter2, int idIter1, int idIter2) :
		UnionTieIterator(iter1, iter2, idIter1, idIter2) {
	init();
}

SymDiffTieIterator::SymDiffTieIterator(const SymDiffTieIterator& rhs) :
		UnionTieIterator(rhs) {
}

SymDiffTieIterator::~SymDiffTieIterator() {
}

SymDiffTieIterator* SymDiffTieIterator::clone() const {
	return new SymDiffTieIterator(*this);
}

void SymDiffTieIterator::init() {
	if (lpIter1->valid() && lpIter2->valid() && isCommon()) {
		next();
	}
}

void SymDiffTieIterator::next() {
	do {
		UnionTieIterator::next();
	} while (lpIter1->valid() && lpIter2->valid() && isCommon());
}

} /* namespace siena */
