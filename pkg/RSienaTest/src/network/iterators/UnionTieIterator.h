/*
 * UnionTieIterator.h
 *
 *  Created on: 25.07.2013
 *      Author: ortmann
 */

#ifndef UNIONTIEITERATOR_H_
#define UNIONTIEITERATOR_H_

#include "CombinedTieIterator.h"

namespace siena {

template<class __iter1, class __iter2>
class UnionTieIterator: public CombinedTieIterator<__iter1, __iter2 > {

public:
	UnionTieIterator(__iter1 iter1, __iter2 iter2) :
			CombinedTieIterator<__iter1, __iter2 >(iter1, iter2) {
	}

	virtual ~UnionTieIterator() {
	}

	inline int actor() const {
		// if both iterators are valid return the actor with the lower value
		if (this->m_iter1.valid() && this->m_iter2.valid()) {
			return std::min(this->m_iter1.actor(), this->m_iter2.actor());
		} else if (this->m_iter1.valid()) {
			return this->m_iter1.actor();
		} else if (this->m_iter2.valid()) {
			return this->m_iter2.actor();
		}
		throw InvalidIteratorException();
	}

	virtual inline void next() {
		// at least one of the iterators has to be valid
		if (!valid()) {
			return;
		}
		// if both iterators are valid move to the next position
		if (this->m_iter1.valid() && this->m_iter2.valid()) {
			if (this->m_iter1.actor() < this->m_iter2.actor()) {
				this->m_iter1.next();
			} else if (this->m_iter1.actor() > this->m_iter2.actor()) {
				this->m_iter2.next();
			} else {
				this->m_iter1.next();
				this->m_iter2.next();
			}
		} else if (this->m_iter1.valid()) {
			this->m_iter1.next();
		} else {
			this->m_iter2.next();
		}
	}

	inline bool valid() const {
		return this->m_iter1.valid() || this->m_iter2.valid();
	}
};

} /* namespace siena */

#endif /* UNIONTIEITERATOR_H_ */
