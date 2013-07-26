/*
 * CombinedTieIterator.h
 *
 *  Created on: 25.07.2013
 *      Author: ortmann
 */

#ifndef COMBINEDTIEITERATOR_H_
#define COMBINEDTIEITERATOR_H_

#include "ITieIterator.h"

namespace siena {

// currently not save ... will be solved once concepts become a standard
template<class __iter1, class __iter2>
class CombinedTieIterator: public ITieIterator {
public:
	virtual ~CombinedTieIterator() {
	}

protected:
	CombinedTieIterator(__iter1 iter1, __iter2 iter2) :
			ITieIterator(), m_iter1(iter1), m_iter2(iter2) {
	}

	inline bool isCommon() {
		return m_iter1.actor() == m_iter2.actor();
	}
	__iter1 m_iter1;
	__iter2 m_iter2;
};

} /* namespace siena */
#endif /* COMBINEDTIEITERATOR_H_ */
