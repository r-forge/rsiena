/*
 * ITieIterator.h
 *
 *  Created on: 16.07.2013
 *      Author: ortmann
 */

#ifndef ITIEITERATOR_H_
#define ITIEITERATOR_H_

#include "../../utils/Utils.h"

namespace siena {

class ITieIterator {
public:
	virtual ~ITieIterator() {
	}
	virtual void next() = 0;
	virtual int actor() const = 0;
	virtual bool valid() const = 0;
	virtual ITieIterator* clone() const = 0;
protected:
	ITieIterator() {
	}

	ITieIterator(const ITieIterator&) {
	}

	ITieIterator& operator=(const ITieIterator&) {
		return *this;
	}
};

} /* namespace siena */
#endif /* ITIEITERATOR_H_ */
