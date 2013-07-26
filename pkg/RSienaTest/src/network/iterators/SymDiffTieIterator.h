/*
 * SymDiffTieIterator.h
 *
 *  Created on: 25.07.2013
 *      Author: ortmann
 */

#ifndef SYMDIFFTIEITERATOR_H_
#define SYMDIFFTIEITERATOR_H_

#include "UnionTieIterator.h"

namespace siena {

template<class __iter1, class __iter2>
class SymDiffTieIterator: public UnionTieIterator<__iter1, __iter2 > {

public:
	SymDiffTieIterator(__iter1 iter1, __iter2 iter2) :
			UnionTieIterator<__iter1, __iter2 >(iter1, iter2) {
		if (this->m_iter1.valid() && this->m_iter2.valid()
				&& this->isCommon()) {
			next();
		}
	}

	~SymDiffTieIterator() {
	}

	inline void next() {
		do {
			UnionTieIterator<__iter1, __iter2 >::next();
		} while (this->m_iter1.valid() && this->m_iter2.valid()
				&& this->isCommon());
	}
};

} /* namespace siena */

#endif /* SYMDIFFTIEITERATOR_H_ */
