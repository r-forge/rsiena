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

class IntersectionTieIterator: public CombinedTieIterator {
public:
	IntersectionTieIterator(const ITieIterator& iter1,
			const ITieIterator& iter2);

	virtual ~IntersectionTieIterator();

	virtual void next();
	virtual int actor() const;
	virtual bool valid() const;
	virtual IntersectionTieIterator * clone() const;
	void skip();

protected:
	IntersectionTieIterator(const IntersectionTieIterator& rhs);
private:
	IntersectionTieIterator& operator=(const IntersectionTieIterator&);

};

} /* namespace siena */

#endif /* INTERSECTIONTIEITERATOR_H_ */
