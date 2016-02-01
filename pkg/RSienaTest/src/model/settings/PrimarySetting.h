/*
 * PrimarySetting.h
 *
 *  Created on: 18.06.2014
 *      Author: ortmann
 */

#ifndef PRIMARYSETTING_H_
#define PRIMARYSETTING_H_

#include "GeneralSetting.h"
#include "../../network/layers/DistanceTwoLayer.h"
namespace siena {

class Network;

class PrimarySetting: public GeneralSetting {
public:
	PrimarySetting();

	virtual ~PrimarySetting();

	/**
	 * @copydoc ASetting::initSetting(Network* const lpNetwork)
	 */
	virtual void initSetting(Network* const lpNetwork);

	/**
	 * @copydoc ASetting::terminateSetting(Network* const lpNetwork)
	 */
	virtual void terminateSetting(Network* const lpNetwork);

	ITieIterator* getSteps();

	int getSize();


protected:

	void initSetting();

	void terminateSetting();

private:

	DistanceTwoLayer rDistTwoLayer;

	Network* lpNetwork;

	ITieIterator* lpiter;

};

} /* namespace siena */
#endif /* PRIMARYSETTING_H_ */
