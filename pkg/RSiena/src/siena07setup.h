/******************************************************************************
 * SIENA: Simulation Investigation for Empirical Network Analysis
 *
 * Web: http://www.stats.ox.ac.uk/~snijders/siena/
 *
 * File: siena07setup.h
 *
 * Description: This file contains prototypes for the Siena data creation
 * funcitons called from R.
 *****************************************************************************/

#ifndef SIENA07SETUP_H_
#define SIENA07SETUP_H_

/**
 *  Creates an array of pointers to Data objects, one for each group
 *  and returns the address of the array to R. Also creates the actor sets
 *  for each group.
 */
SEXP setupData(SEXP OBSERVATIONSLIST, SEXP ACTORSLIST);

/**
 *  Creates all the groups of one mode networks in the data
 *
 */
SEXP OneMode(SEXP RpData, SEXP ONEMODELIST);

/**
 *  Creates all the groups of bipartite networks in the data
 *
 */
SEXP Bipartite(SEXP RpData, SEXP BIPARTITELIST);

/**
 *  Creates all the groups of behavior networks in the data
 */

SEXP Behavior(SEXP RpData, SEXP BEHLIST);

/**
 *  Creates all the groups of constant covariates in the data
 */
SEXP ConstantCovariates(SEXP RpData, SEXP COCOVARLIST);

/**
 *  Creates all the groups of constant covariates in the data
 */
SEXP ChangingCovariates(SEXP RpData, SEXP VARCOVARLIST);

/**
 *  Creates all the groups of constant dyadic covariates in the data
 */
SEXP DyadicCovariates(SEXP RpData, SEXP DYADVARLIST);

/**
 *  Creates all the groups of changing dyadic covariates in the data
 */
SEXP ChangingDyadicCovariates(SEXP RpData, SEXP VARDYADLIST);

/**
 *  Creates all the composition change events in the data
 */
SEXP ExogEvent(SEXP RpData, SEXP EXOGEVENTLIST);

/**
 *  Sets the pairwise constraints for the data
 */
SEXP Constraints(SEXP RpData, SEXP FROMHIGHERLIST, SEXP TOHIGHERLIST,
	SEXP FROMDISJOINTLIST, SEXP TODISJOINTLIST,
	SEXP FROMATLEASTONELIST, SEXP TOATLEASTONELIST);

/**
 *  creates the requested basic effects
 */
SEXP effects(SEXP RpData, SEXP EFFECTSLIST);

/**
 *  creates the requested interaction effects
 */
SEXP interactionEffects(SEXP RpData, SEXP RpModel, SEXP EFFECTSLIST);

/**
 *  removes the objects created for the data.
 */
SEXP deleteData(SEXP RpData);

/**
 *  removes the model object.
 */
SEXP deleteModel(SEXP RpModel);

/**
 *  sets up the model options of MAXDEGREE, CONDITIONAL
 */
SEXP setupModelOptions(SEXP DATAPTR, SEXP MODELPTR, SEXP MAXDEGREE,
	SEXP CONDVAR, SEXP CONDTARGETS, SEXP PROFILEDATA, SEXP PARALLELRUN);

/**
 *  Gets target values relative to the input data
 */
SEXP getTargets(SEXP DATAPTR, SEXP MODELPTR, SEXP EFFECTSLIST);

#endif /*SIENA07SETUP_H_*/