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

class CombinedTieIterator: public ITieIterator {
public:
	virtual ~CombinedTieIterator();

protected:
	CombinedTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);
	CombinedTieIterator(const CombinedTieIterator& rhs);

	inline bool isCommon() const {
		return lpIter1->actor() == lpIter2->actor();
	}

protected:
	ITieIterator* const lpIter1;
	ITieIterator* const lpIter2;
private:
	CombinedTieIterator& operator=(const CombinedTieIterator&);
};

} /* namespace siena */
#endif /* COMBINEDTIEITERATOR_H_ */
