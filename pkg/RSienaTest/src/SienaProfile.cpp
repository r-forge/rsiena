//============================================================================
// Name        : SienaProfile.cpp
// Author      : Ruth Ripley
// Version     :
// Copyright   :
// This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, a copy is available at
//  http://www.r-project.org/Licenses///
//============================================================================

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <stdexcept>
#include <cstdlib>
#include <sys/time.h>
#include "network/OneModeNetwork.h"
#include "network/TieIterator.h"
#include "network/NetworkUtils.h"
#include "data/Data.h"
#include "data/ActorSet.h"
#include "data/LongitudinalData.h"
#include "data/NetworkLongitudinalData.h"
#include "data/OneModeNetworkLongitudinalData.h"
#include "data/BehaviorLongitudinalData.h"
#include "data/ChangingDyadicCovariate.h"
#include "data/ConstantDyadicCovariate.h"
#include "data/ChangingCovariate.h"
#include "data/ConstantCovariate.h"
#include "data/ExogenousEvent.h"

#include "model/Model.h"
#include "model/EpochSimulation.h"
#include "model/State.h"
#include "model/StatisticCalculator.h"
#define MATHLIB_STANDALONE
#include <Rmath.h>

using namespace siena;
using namespace std;

double grandSum;
double scoreSum;

vector<Data *> * readInData( )
{
	ifstream myfile ("data.txt");

	vector<Data *> *pGroupData = new vector <Data *>;

	int nGroups = 1;
	int observations;
	int nActors;
	myfile >> observations;
	myfile >> nActors;

	const ActorSet *pActors;
	OneModeNetworkLongitudinalData *  pOneModeNetworkLongitudinalData;


	for (int group = 0; group < nGroups; group++)
	{

		pGroupData->push_back(new Data(observations));

		int nNodeSets = 1;
		for (int nodeSet = 0; nodeSet < nNodeSets; nodeSet++)
		{

			pActors =	(*pGroupData)[group]->
				createActorSet("Actors", nActors);
		}
	}

	string networkName;
	string variableType;
	Data * pData = (*pGroupData)[0];
	if (myfile.is_open())
	{
		myfile >> networkName >> variableType;
		cout << networkName << endl;
		pOneModeNetworkLongitudinalData =
			pData->createOneModeNetworkData(networkName, pActors);
		pOneModeNetworkLongitudinalData->symmetric(FALSE);
		for (int period = 0; period < (observations - 1); period++)
		{
			pOneModeNetworkLongitudinalData->upOnly(period,
					FALSE);
			pOneModeNetworkLongitudinalData->downOnly(period,
					FALSE);
		}
		for (int period = 0; period < observations; period++)
		{
			int tieCount = 0;
			myfile >> tieCount;
			//	cout << tieCount << endl;
			int ego;
			int alter;
			double value;
			for (int tie = 0; tie < tieCount; tie++)
			{
				myfile >> ego >> alter >> value;
				//				cout << ego << " "<< alter << " " << value << endl;
				pOneModeNetworkLongitudinalData->tieValue(ego,
					alter,
					period,
					value);
			}
			myfile >> tieCount;
			//	cout << "MISSING " <<tieCount << endl;
			for (int tie = 0; tie < tieCount; tie++)
			{
				//	cout << "here " << endl;
				myfile >> ego >> alter >> value;
				//			cout << ego << " "<< alter << " " << value << endl;
				pOneModeNetworkLongitudinalData->missing(ego,
					alter,
					period,
					value);
			}

			myfile >> tieCount;
			//	cout <<"structr tie " << tieCount << endl;
			for (int tie = 0; tie < tieCount; tie++)
			{
				//	cout << "here " << endl;
				myfile >> ego >> alter >> value;
				//				cout << ego << " "<< alter << " " << value << endl;
				pOneModeNetworkLongitudinalData->structural(ego,
					alter,
					period,
					value);
			}
			if (period < observations)
			{
				int uponly;
				myfile >> uponly;
				//	cout << "up " << uponly << endl;
				pOneModeNetworkLongitudinalData->upOnly(period,
					uponly);
				int downonly;
				myfile >> downonly;
				//	cout << "down " << downonly << endl;
				pOneModeNetworkLongitudinalData->downOnly(period,
					downonly);
			}
		}

		// Once all network data has been read, calculate some
		// statistical properties of that data.

		pOneModeNetworkLongitudinalData->calculateProperties();

		// may be behavior variables etc
		string variableName;
		int n;
		int actor;
		int alter;
		int value;
		double realvalue;
		double simvalue;
		double range;
		double mean;
		myfile >> variableName >> variableType;
		cout << variableName << " " << variableType << endl;
		while (!myfile.eof())
		{
			if (variableType =="behavior")
			{
				BehaviorLongitudinalData * pBehaviorData =
					pData->createBehaviorData(variableName, pActors);
				myfile >> n;
				//	cout << n << endl;
				for (int period = 0; period < observations; period++)
				{
					for (int i = 0; i < n; i++)
					{
						myfile >> actor >> realvalue;
						//			cout << actor << " " << realvalue << endl;
						pBehaviorData->value(period, actor, realvalue);
					}
					for (int i = 0; i < n; i++)
					{
						myfile >> actor >> value;
						//			cout << actor << " " << value << endl;
						pBehaviorData->missing(period, actor, value);
					}
					myfile >> simvalue;
					//		cout << simvalue << endl;
					pBehaviorData->similarityMean(simvalue);
					if (period < observations)
					{
						int uponly;
						myfile >> uponly;
						//		cout << uponly << endl;
						pBehaviorData->upOnly(period, uponly);
						int downonly;
						myfile >> downonly;
						//		cout << downonly << endl;
						pBehaviorData->downOnly(period, downonly);
					}
					// Now that the values are set,
					// calculate some important statistics
					pBehaviorData->calculateProperties();
				}
			}
			else if (variableType == "constantcovariate")
			{
				ConstantCovariate * pConstantCovariate =
					pData->createConstantCovariate(variableName,
						pActors);
				myfile >> n;
				//	cout << n << endl;

				for (int i = 0; i < n; i++)
				{
					myfile >> actor >> realvalue;
					//				cout << actor << " " << realvalue << endl;
					pConstantCovariate->value(actor, realvalue);
				}
				for (int i = 0; i < n; i++)
				{
					myfile >> actor >> value;
					//	cout << actor << " " << value << endl;
					pConstantCovariate->missing(actor, value);
				}
				myfile >> simvalue;
				//	cout << simvalue << endl;
				pConstantCovariate->similarityMean(simvalue);
				myfile >> range;
				//	cout << range << endl;
				pConstantCovariate->range(range);
			}
			else if (variableType == "changingcovariate")
			{
				ChangingCovariate * pChangingCovariate =
					pData->createChangingCovariate(variableName,
						pActors);
				myfile >> n;
				//	cout << n << endl;
				for (int period = 0; period < observations - 1; period++)
				{

					for (int i = 0; i < n; i++)
					{
						myfile >> actor >> realvalue;
						//		cout << actor << " " << realvalue << endl;
						pChangingCovariate->value(actor, realvalue, period);
					}
					for (int i = 0; i < n; i++)
					{
						myfile >> actor >> value;
						//		cout << actor << " " << value << endl;
						pChangingCovariate->missing(actor, value, period);
					}
				}
				myfile >> simvalue;
				//	cout << simvalue << endl;
				pChangingCovariate->similarityMean(simvalue);
				myfile >> range;
				//	cout << range << endl;
				pChangingCovariate->range(range);

			}
			else if (variableType == "constantdyadiccovariate")
			{
				ConstantDyadicCovariate * pCovariate =
					pData->createConstantDyadicCovariate(variableName,
						pActors, pActors);
				myfile >> n;
				//	cout << n << endl;
				for (int i = 0; i < n; i++)
				{
					for (int j = 0; j < n; j++)
					{
						myfile >> actor >> alter >> value;
						//cout << actor << " " << alter << " " << value << endl;
						pCovariate->value(actor, alter, value);
						myfile >> actor >> alter >> value;
						//cout << actor << " " << alter << " " << value << endl;
						pCovariate->missing(actor, alter, value);
					}
				}
				myfile >> mean;
				//	cout << mean << endl;
				pCovariate->mean(mean);
			}
			else if (variableType == "changingdyadiccovariate")
			{
				ChangingDyadicCovariate * pCovariate =
					pData->createChangingDyadicCovariate(variableName,
						pActors, pActors);
				myfile >> n;
				//	cout << n << endl;

				for (int period = 0; period < observations - 1;
					 period++)
				{
					for (int i = 0; i < n; i++)
					{
						for (int j = 0; j < n; j++)
						{
							myfile >> actor >> alter >> value;
							//cout << actor << " " << alter << " " <<
							// value << endl;
							pCovariate->value(actor, alter,	period,	value);
							myfile >> actor >> alter >> value;
							//cout << actor << " " << alter << " " << value
							//<< endl;
							pCovariate->missing(actor, alter, period, value);
						}
					}
				}
				myfile >> mean;
				//	cout << mean << endl;
				pCovariate->mean(mean);
			}
			else
			{
						cout << " write another module" << endl;
			}
			myfile >> variableName >> variableType;
			cout << variableName << " " << variableType << endl;
		}
	}
	else
	{
		cout << "not open\n";
	}
	return pGroupData;
}
Model *  readInEffects(vector<Data *> *pGroupData)
{
	int nGroups = pGroupData->size();
	int totObservations = 0;
	for (int group = 0; group < nGroups; group++)
		totObservations += (*pGroupData)[group]->observationCount() - 1;

	Model * pModel = new Model();

	string variableName;
	string effectName;
	string effectType;
	double parameter;
	double internalEffectParameter;
	string interactionName1;
	string interactionName2;
	string rateType;
	string netType;
	int period;
	int group;
	string mystring;
	istringstream ss;
	ifstream effectsfile ("effects.txt");
	if (effectsfile.is_open())
	{
		while (getline(effectsfile, mystring))
		{
			ss.clear();
			ss.str(mystring);
			ss >> variableName >> effectName >> effectType>> parameter
			   >>internalEffectParameter >>
				interactionName1 >>interactionName2 >>
				rateType >> netType >> group ;
			cout << variableName <<" " << effectName <<" " <<
				effectType<< " "<< parameter << " "
				 <<internalEffectParameter <<" " << interactionName1<<" "<<
				interactionName2 << " " << rateType << " " <<
				netType <<" " << group << "group" << endl;
			if (interactionName1 == "NA")
			{
				interactionName1 = "";
			}
			if (interactionName2 == "NA")
			{
				interactionName2 = "";
			}
			if (rateType == "NA")
			{
				rateType = "";
			}
			if (effectName == "Rate")
			{
				Data *pData = (*pGroupData)[group-1];
				ss >> period;
				cout << "period " << period << endl;
				if (netType == "oneMode")
				{
					//	cout << "here\n" << (void *) pGroupData;
					NetworkLongitudinalData * pNetwork =
						pData->pNetworkData(variableName);
					//	cout << "here1\n";
					pModel->basicRateParameter(pNetwork, period-1, parameter);
					//	cout << "here2\n";
				}
				else if (netType == "Behavior")
				{
					BehaviorLongitudinalData * pNetwork =
						pData->pBehaviorData(variableName);
					pModel->basicRateParameter(pNetwork, period-1, parameter);
				}
			}
			else
			{
				pModel->addEffect(variableName, effectName, effectType,
					parameter, internalEffectParameter,
					interactionName1, interactionName2, rateType);
			}
		}
		effectsfile.close();
	}
	else
	{
		cout << "effects file not open";
	}
	return pModel;
}
void FitModel(vector<Data *> * pGroupData, Model * pModel, int deriv)
{
	/* create a simulation and return the observed statistics and scores */

	int nGroups = pGroupData->size();

	int totObservations = 0;
	for (int group = 0; group < nGroups; group++)
		totObservations += (*pGroupData)[group]->observationCount() - 1;

	int periodFromStart = 0;
    /* group loop here */
	for (int group = 0; group < nGroups; group++)
	{

		/* find out how many periods in this Data object */
		Data * pData = (*pGroupData)[group];
		int observations = pData->observationCount();

		/* create my epochsimulation object */
		EpochSimulation * pEpochSimulation  = new
		EpochSimulation(pData, pModel);

		for (int period = 0; period < observations - 1; period++)
        {
			periodFromStart++;

			/* run the epoch simulation for this period */
			pEpochSimulation->runEpoch(period);
			State State(pEpochSimulation);
			StatisticCalculator Calculator(pData, pModel, &State,
				period);

			for (unsigned i = 0;
				i < pData->rDependentVariableData().size();
				i++)
			{
				LongitudinalData * pVariableData =
					pData->rDependentVariableData()[i];
				string name = pVariableData->name();

				for (unsigned j = 0;
					j < pModel->rEvaluationEffects(name).size();
					j++)
				{
					EffectInfo * pInfo = pModel->rEvaluationEffects(name)[j];
					grandSum += Calculator.statistic(pInfo);
					scoreSum += pEpochSimulation->score(pInfo);
				}

				for (unsigned j = 0;
					j < pModel->rEndowmentEffects(name).size();
					j++)
				{
					EffectInfo * pInfo = pModel->rEndowmentEffects(name)[j];
					grandSum += Calculator.statistic(pInfo);
					scoreSum += pEpochSimulation->score(pInfo);
				}
			}
        } /* end of period */
			delete pEpochSimulation;
	}	 /* end of group */

    return ;
}



int main(int argc, char** argv) {

	int iterations = atoi(argv[1]);
	int deriv = atoi(argv[2]);

	cout << "deriv " << deriv << " number of iterations " << iterations << endl;
	time_t rawtime, endtime;
	time ( &rawtime );
	cout << "The current local time is: " << ctime (&rawtime) << endl ;

	vector<Data *> * pGroupData = readInData();
	cout << "got here"<<endl;
	Model * pModel = readInEffects(pGroupData);
	pModel->needScores(deriv);
	cout << "model need scores "<< pModel->needScores() << " deriv " << deriv << endl;
	time ( &rawtime );
	cout << "The current local time is: " << ctime (&rawtime) << endl ;
	grandSum = 0;
	scoreSum = 0;
	for (int i=0; i < iterations; i++)
	{
		//	cout << i<< endl;
		FitModel(pGroupData, pModel, deriv);
	}
	cout << "Checksum: " << grandSum << " score "<< scoreSum << endl;
	time ( &endtime );
	double diff = endtime - rawtime;
	cout << "The current local time is: " << ctime (&endtime) << endl ;
	cout << "time taken per it: " << diff/iterations << endl ;
	cout << "OK" << endl;
	return EXIT_SUCCESS;
}
