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

class SymDiffTieIterator: public UnionTieIterator {
public:
	SymDiffTieIterator(const ITieIterator& iter1, const ITieIterator& iter2);
	SymDiffTieIterator(const ITieIterator& iter1, const ITieIterator& iter2,
			int idIter1, int idIter2);
	virtual ~SymDiffTieIterator();
	virtual SymDiffTieIterator* clone() const;
	virtual void next();
protected:
	SymDiffTieIterator(const SymDiffTieIterator& rhs);
private:
	void init();
};

} /* namespace siena */

#endif /* SYMDIFFTIEITERATOR_H_ */
