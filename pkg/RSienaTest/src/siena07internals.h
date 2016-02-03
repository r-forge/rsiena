/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07internals.h
 *
 * Description: This file contains prototypes for the internal routines
 * used to setup the data in C and create initial ML chains.
 *****************************************************************************/

#ifndef SIENA07INTERNALS_H_
#define SIENA07INTERNALS_H_

#include <Rinternals.h>

#include <vector>

namespace siena
{

	class Data;
	class Model;
	class StatisticCalculator;
	class EpochSimulation;
	class MLSimulation;
	class NetworkLongitudinalData;
	class OneModeNetworkLongitudinalData;
	class BehaviorLongitudinalData;
	class ConstantCovariate;
	class ChangingCovariate;
	class ConstantDyadicCovariate;
	class ChangingDyadicCovariate;

}

/**
 * Matches column names with indices in the effects object.
 */
void getColNos(SEXP Names, int * netTypeCol, int * nameCol, int * effectCol,
	int * parmCol, int * int1Col, int * int2Col, int * initValCol,
	int * typeCol, int * groupCol, int * periodCol, int * pointerCol,
	int * rateTypeCol, int * intptr1Col, int * intptr2Col, int * intptr3Col,
	int * settingCol);

/**
 *  updates the parameter values for each of the effects.
 */
void updateParameters(SEXP EFFECTSLIST, SEXP THETA, std::vector<siena::Data *>
		*pGroupData, siena::Model *pModel);

/**
 * Create one observation for a one mode Network: ties, missing, structural
 *
 */
void setupOneModeNetwork(SEXP ONEMODE, siena::OneModeNetworkLongitudinalData
		*pNetworkData, int observation);

/**
 * Create all observations for a one mode Network
 *
 */
void setupOneModeObservations(SEXP ONEMODES,
	siena::OneModeNetworkLongitudinalData *pOneModeNetworkLongitudinalData);

/**
 * Create one group of one mode Networks
 *
 */
void setupOneModeGroup(SEXP ONEMODEGROUP, siena::Data * pData);

/**
 * Create one observation for a bipartite Network: ties, missing, structural
 *
 */
void setupBipartiteNetwork(SEXP BIPARTITE,
	siena::NetworkLongitudinalData * pNetworkData,
	int observation);

/**
 * Create all observations for a bipartite Network
 *
 */
void setupBipartiteObservations(SEXP BIPARTITES,
	siena::NetworkLongitudinalData *
	pNetworkLongitudinalData);

/**
 * Create one group of bipartite Networks
 *
 */
void setupBipartiteGroup(SEXP BIPARTITEGROUP, siena::Data *pData);

/**
 * Create all observations for a behavior Network
 *
 */
void setupBehavior(SEXP BEHAVIOR, siena::BehaviorLongitudinalData *pBehaviorData);


/**
 * Create one group of Behavior Networks
 *
 */
void setupBehaviorGroup(SEXP BEHGROUP, siena::Data *pData);

/**
 * Create a constant covariate
 *
 */
void setupConstantCovariate(SEXP COCOVAR, siena::ConstantCovariate
		*pConstantCovariate);

/**
 * Create one group of constant covariates
 *
 */
void setupConstantCovariateGroup(SEXP COCOVARGROUP, siena::Data *pData);

/**
 * Create all observations for a changing covariate
 *
 */
void setupChangingCovariate(SEXP VARCOVAR,
		siena::ChangingCovariate * pChangingCovariate);


/**
 * Create one group of changing covariates
 *
 */
void setupChangingCovariateGroup(SEXP VARCOVARGROUP, siena::Data *pData);

/**
 * Create a constant dyadic covariate
 *
 */
void setupDyadicCovariate(SEXP DYADVAR,
		siena::ConstantDyadicCovariate * pConstantDyadicCovariate);

/**
 * Create one group of constant dyadic covariates
 *
 */
void setupDyadicCovariateGroup(SEXP DYADVARGROUP, siena::Data *pData);

/**
 * Unpack one set of values for a changing dyadic covariate
 *
 */
void unpackChangingDyadicPeriod(SEXP VARDYADVALS,
		siena::ChangingDyadicCovariate *pChangingDyadicCovariate, int period);

/**
 * Create all observations for a changing dyadic covariate
 *
 */
void setupChangingDyadicObservations(SEXP VARDYAD,
		siena::ChangingDyadicCovariate *pChangingDyadicCovariate);

/**
 * Create one group of changing dyadic covariates
 *
 */
void setupChangingDyadicCovariateGroup(SEXP VARDYADGROUP, siena::Data * pData);

/**
 * Create the exogenous composition change events for one actor set within
 * one group.
 */
void setupExogenousEventSet(SEXP EXOGEVENTSET, siena::Data *pData);

/**
 * Create one group of exogenous composition change events
 *
 */
void setupExogenousEventGroup(SEXP EXOGEVENTGROUP, siena::Data *pData);

/**
 *  Creates all the basic effects for one network
 */
SEXP createEffects(SEXP EFFECTS, siena::Model *pModel,
		std::vector<siena::Data *> *pGroupData,
		const char *networkName, int effectCol,
		int parmCol, int int1Col, int int2Col,
		int initValCol, int typeCol, int groupCol,
		int periodCol, int rateTypeCol,
		int netTypeCol, int settingCol);

/**
 *  Creates all the interaction effects for one network
 */
SEXP createInteractionEffects(SEXP EFFECTS, siena::Model *pModel,
		const char *networkName, int effectCol, int initValCol,
		int typeCol, int intptr1Col, int intptr2Col, int intptr3Col);

/**
 *  Retrieves the contributions to all possible tie flips or behavior changes for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getChangeContributionStatistics(SEXP EFFECTSLIST,
		const siena::StatisticCalculator *pCalculator,
		std::vector<std::vector<double * > > *rChangeContributions);

/**
 *  Retrieves the statistics of individual actors for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getActorStatistics(SEXP EFFECTSLIST,
	const siena::StatisticCalculator * pCalculator, std::vector<double *> *rActorStatistics);

/**
 *  Retrieves the values of the statistics and scores for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Not used in maximum likelihood.
 */
void getStatistics(SEXP EFFECTSLIST,
		const siena::StatisticCalculator *pCalculator,
		int period, int group, const siena::Data *pData,
		const siena::EpochSimulation *pEpochSimulation,
		std::vector<double> * rfra, std::vector<double> *rscore);

/**
 *  retrieves the values of the scores and derivatives for each of the effects,
 *  for one period. The call will relate to one group only, although all effects
 *  are the same apart from the basic rates. Only used in maximum likelihood.
 */
void getScores(SEXP EFFECTSLIST, int period, int group,
		const siena::MLSimulation *pMLSimulation,
		std::vector<double> *rderiv, std::vector<double> *rscore);

#endif /*SIENA07INTERNALS_H_*/
