#/******************************************************************************
# * SIENA: Simulation Investigation for Empirical Network Analysis
# *
# * Web: http://www.stats.ox.ac.uk/~snijders/siena
# *
# * File: tests/slowtest.R
# *
# * Description: This module tests the examples on the help pages,
# * including those marked dontrun.
# *
# * Not run automatically, as excluded from the tar ball.
# *****************************************************************************/
## R CMD BATCH --no-save "--args RSienaDir='~/rforge/pkg/RSienaTest'
## lib='library(RSienaTest)'" ~/rforge/pkg/RSienaTest/tests/slowtest.R
## (without the carriage return)
## First argument is where to look for Rd files
## Second one is a command to be executed before running the code.
## At least this needs to load either RSiena or RSienaTest, but could do more.
## Note that there are functions in RSienaTest which do not exist in RSiena.
## Can be run inside R if you prefer.
xx <- commandArgs(trailing=TRUE)
eval(parse(text=xx))
RSienaTest:::RSienaSlowTest(RSienaDir, lib)
