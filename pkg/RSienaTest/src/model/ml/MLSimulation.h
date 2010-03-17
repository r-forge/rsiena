#ifndef MLSIMULATION_H_
#define MLSIMULATION_H_

#include "model/EpochSimulation.h"

namespace siena
{

// ----------------------------------------------------------------------------
// Section: Forward declarations
// ----------------------------------------------------------------------------

class Chain;
class MiniStep;


// ----------------------------------------------------------------------------
// Section: Class definition
// ----------------------------------------------------------------------------

class MLSimulation: public EpochSimulation
{
public:
	MLSimulation(Data * pData, Model * pModel);
	virtual ~MLSimulation();

	void connect(int period);
    void updateProbabilities(Chain * pChain,
    	MiniStep * pFirstMiniStep,
    	MiniStep * pLastMiniStep);
    void executeMiniSteps(MiniStep * pFirstMiniStep, MiniStep * pLastMiniStep);

    // Metropolis-Hastings steps

	bool insertDiagonalMiniStep();
	bool cancelDiagonalMiniStep();
	bool permute(int c0);
	double proposalProbability() const;

	void simpleRates(bool flag);
	bool simpleRates() const;

	void insertDiagonalProbability(double probability);
	double insertDiagonalProbability() const;

	void cancelDiagonalProbability(double probability);
	double cancelDiagonalProbability() const;

private:
	void setStateBefore(MiniStep * pMiniStep);
	void resetVariables();

	Chain * lpChain;
	bool lsimpleRates;
	double linsertDiagonalProbability;
	double lcancelDiagonalProbability;
	double lproposalProbability;
};

}

#endif /* MLSIMULATION_H_ */
