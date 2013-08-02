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

class UnionTieIterator: public CombinedTieIterator {
public:
	UnionTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);
	UnionTieIterator(const ITieIterator& iter1, const ITieIterator& iter2,
			int idIter1, int idIter2);
	virtual ~UnionTieIterator();

	virtual void next();
	int actor() const;
	bool valid() const;
	virtual UnionTieIterator * clone() const;

	int getInactiveIterID();
	int getActiveIterID();
	bool isCommonNeighbor() const;

protected:
	UnionTieIterator(const UnionTieIterator& rhs);

private:
	int lIdIter1;
	int lIdIter2;

	void init(int idIter1, int idIter2);
	UnionTieIterator& operator=(const UnionTieIterator&);
};

} /* namespace siena */

#endif /* UNIONTIEITERATOR_H_ */
