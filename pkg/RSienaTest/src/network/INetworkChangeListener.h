/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: NetworkChangeListener.h
 *
 * Description: This module defines the interface INetworkChangeListener.
 * Any class implementing this interface can be added to a network and gets
 * informed once an edge is introduced/withdrawn from the network.
 *****************************************************************************/

#ifndef INETWORKCHANGELISTENER_H_
#define INETWORKCHANGELISTENER_H_

#include "Network.h"

namespace siena {

class INetworkChangeListener {
public:
	virtual ~INetworkChangeListener() {
	}
	virtual void onTieIntroductionEvent(const Network& rNetwork, const int ego,
			const int alter) = 0;
	virtual void onTieWithdrawalEvent(const Network& rNetwork, const int ego,
			const int alter) = 0;
	virtual void onNetworkClearEvent(const Network& rNetwork) = 0;
protected:
	INetworkChangeListener() {
	}
private:
	// disable copy constructor and copy assignment
	INetworkChangeListener(const INetworkChangeListener& rhs);
	INetworkChangeListener& operator=(const INetworkChangeListener& rhs);
};

}
#endif /* INETWORKCHANGELISTENER_H_ */
