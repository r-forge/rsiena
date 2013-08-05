/*
 * NetworkLayer.h
 *
 *  Created on: 05.08.2013
 *      Author: ortmann
 */

#ifndef NETWORKLAYER_H_
#define NETWORKLAYER_H_

#include "../INetworkChangeListener.h"

namespace siena {
class NetworkLayer: public INetworkChangeListener {
public:
	virtual ~NetworkLayer() {
	}
protected:
	virtual void initializeOneMode(const Network& rNetwork) = 0;
	NetworkLayer() :
			INetworkChangeListener() {
	}
private:
	// disable copy constructor and copy assignment
	NetworkLayer& operator=(const NetworkLayer& rhs);
	NetworkLayer(const NetworkLayer& rhs);
};

} /* namespace siena */

#endif /* NETWORKLAYER_H_ */
