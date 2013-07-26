/*
 * IntersectionTieIterator.h
 *
 *  Created on: 25.07.2013
 *      Author: ortmann
 */

#ifndef INTERSECTIONTIEITERATOR_H_
#define INTERSECTIONTIEITERATOR_H_

#include "CombinedTieIterator.h"

namespace siena {

template<class __iter1, class __iter2>
class IntersectionTieIterator: public CombinedTieIterator<__iter1, __iter2 > {

public:
	IntersectionTieIterator(__iter1 iter1, __iter2 iter2) :
			CombinedTieIterator<__iter1, __iter2 >(iter1, iter2) {
		// move to the first common position
		if (valid() && !isCommon()) {
			skip();
		}
	}

	~IntersectionTieIterator() {
	}

	inline int actor() const {
		// we only have to check whether iter2 is valid cause
		// iter1 will throw anyway
		if (this->m_iter2.valid()) {
			return this->m_iter1.actor();
		}
		throw InvalidIteratorException();
	}

	inline void next() {
		this->m_iter1.next();
		this->m_iter2.next();
		skip();
	}

	inline void skip() {
		while (valid() && !this->isCommon()) {
			while (this->m_iter1.valid()
					&& this->m_iter1.actor() < this->m_iter2.actor()) {
				this->m_iter1.next();
			}
			if (!this->m_iter1.valid()) {
				return;
			}
			while (this->m_iter2.valid()
					&& this->m_iter2.actor() < this->m_iter1.actor()) {
				this->m_iter2.next();
			}
		}
	}

	inline bool valid() const {
		return this->m_iter1.valid() && this->m_iter2.valid();
	}

};

//typedef IntersectionTieIterator<IncidentTieIterator, IncidentTieIterator> CommonNeighborIterator;

} /* namespace siena */

#endif /* INTERSECTIONTIEITERATOR_H_ */
