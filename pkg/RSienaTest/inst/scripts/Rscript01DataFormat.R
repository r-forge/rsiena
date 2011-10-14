###################################################################################
###
### ---- RscriptDataFormat.R: a script for the introduction to RSiena -------------
###
###                               version: April 11, 2011
###################################################################################
#
# RscriptDataFormat.R is followed by
# RScriptSNADescriptives.R, code for descriptive analysis of the data, and
# RscriptSienaVariableFormat.R, which formats data and specifies the model, and
# RscriptSienaRunModel.R, which runs the model and estimates parameters
# RscriptSienaBehaviour.R, which illustrates an example of analysing the
# coevolution of networks and behaviour
#
# The entire model fitting is summarised at the end of RscriptSienaRunModel.R
# (without comments)
#
# This is an R script for getting started with RSiena, written by
# Robin Gauthier, Tom Snijders, Ruth Ripley, Johan Koskinen, and
# Paulina Preciado, with some examples borrowed from Christian Steglich.
# Lines starting with # are not processed by R but treated as comments.
# The script has a lot of explanation of R possibilities that will be
# familiar for readers well acquainted with R, and can be skipped by them.
#
# A terrifically easily accessible online introduction to R is given at
#
#        http://www.statmethods.net/
#
# R is case sensitive.
# The left-arrow "<-" is very frequently used: it denotes an assignment,
# "a <- b" meaning that object a gets the value b.
# Often b is a complicated expression that has to be evaluated by R.

# Help within R can be called by typing a question mark and the name of the
# function you need help with. For example ?library loading will bring up a
# file titled "loading and listing of packages".
# Comments are made at the end of commands after #,
# or in lines staring with # telling R to ignore everything beyond it.
# This session will be using s50 data which are supposed to be
# present in the working directory.
# Note that any command in R is called a function;
# in general the command syntax for calling R's functions is function(x) where
# function is a saved function and x the name of the object to be operated on.
# For new R users:
# note that there is a lot of documentation available at
# http://cran.xl-mirror.nl/other-docs.html
# including some short introductions, handy reference cards,
# and introductions in a lot of languages besides English.

#################### - CALLING THE DATA AND PRELIMINARY MANIPULATIONS - ###########

# The library command loads the packages needed during the session.

        library(RSiena)
        library(snow) # (these four additional libraries will be loaded
        library(network)# automatically if required)
        library(rlecuyer)

# You need to have INSTALLED all of them

	?install.packages

# Or click on the tab "Packages", "Instal package(s)", then select a CRAN mirror
# (e.g. Bristol if you are in the UK) and finally select from the list
# the package you wish to install.

# Where are you?

	getwd()

# By something like
	# setwd('C:/SienaTest')
# you can set the directory but note the quotes and forward slash.
# It is also possible to set the directory using the menus if you have them.
# On a windows machine, you can predetermine the working directory
# in the <Properties> of the shortcut to R;
# these are accessible by right-clicking the shortcut.

# If you want to set the working directory to for example
	# "C:\Documents and Settings\johan\My Documents\RSiena_course"
# simply copy and paste from windows explorer or type
#     setwd('C:/Documents and Settings/johan/My Documents/RSiena_course')
# or
#     setwd('C:\\Documents and Settings\\johan\\My Documents\\RSiena_course')
#
# but note that "\" has to be changed; both '/' and '\\' work!!!

# What is there?

         list.files()

# What is available in RSiena?

         ?RSiena

# (these are .htlm HELP PAGES)

# Or, for a listing of all accessible functions in RSiena:

         library(help=RSiena)

# which is very useful

# Where is the manual?

         RShowDoc("RSiena_Manual", package="RSiena")

# (Note, however, that it is possible that the Siena website
# at http://www.stats.ox.ac.uk/~snijders/siena/ contains a more recent version.)

# Each data is named (for example below we name it friend.data.w1)
# so that we can call it as an object within R.
# If you read an object straight into R, it will treat it as a
# dataset, or in R terminology a "data frame".
# Here this is not what we want, therefore on reading
# we will immediately convert it to a matrix.
# R will read in many data formats, these are saved as .dat files, the command
# to read them is read.table.
# If we wished to read a .csv file we would have
# used the read.csv command.
# The pathnames have forward slashes, or double backslashes
# if single backslashes are used, one of the error messages will be:
#   1: '\R' is an unrecognized escape in a character string

#-----------------------------------------------------------------------------

# Quick start(Data assignment).
# Please make sure the data is in your working directory

        friend.data.w1 <- as.matrix(read.table("s50-network1.dat"))
        friend.data.w2 <- as.matrix(read.table("s50-network2.dat"))
        friend.data.w3 <- as.matrix(read.table("s50-network3.dat"))
        drink <- as.matrix(read.table("s50-alcohol.dat"))
        smoke <- as.matrix(read.table("s50-smoke.dat"))

# You need the data to be in your working directory.
# The data set is on ther Siena website ("Datasets" tab) and must be
# unzipped in your working directory.

# Explanation of data structures and formats below

################# - DIFFERENT DATA FORMATS - #######################################
###
### The assignments here involve reading in data to a "data frame"
### 		data <- read.table("data.dat")
### reads your text file into a "data frame"; check the class of the object "data"
###        >  class(data)
###        [1] "data.frame"
### If your data is not a ".dat" file, alternative import methods are
### ----- ".csv" ------------------------------------------------------------------
###        data <- read.table("c:/data.csv", header=TRUE,
###                                   sep=",", row.names="iden")
### ---- ".csv" -------------------------------------------------------------------
###        data <- read.csv("data.csv",header=T)
### ---- ".spss" ------------------------------------------------------------------
###        library(foreign)
###        data <- read.spss("c:/data.spss")
### ---- ".dta" -------------------------------------------------------------------
###        library(foreign)
###        data <- read.dta("c:/data.dta")
###
####################################################################################

################# - FROM DATA FRAME TO MATRIX - ####################################
###
### ---- Data frame ----------------------------------------------------------------
### The data frame is like a spreadsheet of cases by variables,
### where the variables are the columns, and these have the names
###        names( data )
### a spreadsheet view of a data frame is given by the fix() command
###        fix( data )
###
### All changes you make in this spreadsheet will be saved automatically,
### so beware.
###
### As an example create two vectors:
	      height <- c( 120, 150, 180, 200, 190, 177, 170, 125, 141, 157 )
      	weight <- c( 11, 14, 17, 18, 17, 18, 11, 12, 10, 15 )
### collect these in a data frame
	      data <- data.frame( height, weight )
### and look at the results
		    data
### The columns of a data frame may be extracted using a "$" sign and their names.
### For example:
		names( data )
		data$height
### Or by "[]" and column number, e.g.
		data[1]
  	data[  , 1]
  	data[ 1,  ]

### ---- Matrix -------------------------------------------------------------------
### A "matrix" is a numeric object ordered like a matrix with dimensions
### ( dim() ) given by the number of rows and columns, e.g.
		dim( friend.data.w1 )

### If you wish to play around with a copy of the matrix, e.g. having the name "data",
### you can make the copy by the command
		data <- friend.data.w1
### the element of a matrix called data that is in row 2 and column 3 is given by
		data[ 2, 3 ]
### the first three rows of a matrix called data are given by
		data[ 1:3, ]
### columns 2, 5, and 6 are given by
		data[ , c( 2, 5, 6) ]

### ---- Converting data frame to matrix -------------------------------------------
###
### Most classes can be converted to other classes through "as.typeofclass ()", e.g.
	      data <- as.matrix( data )
### converts a data frame "data" to a matrix "data"

####################################################################################


################## - EXAMPLE FOR ARC LIST - ########################################
###
### If data (e.g. "s50e.dat") is in arclist format, with columns
### sender id, receiver id, value of tie, wave.
		ArcList <- read.table( "s50e.dat", header=FALSE ) # creates data frame
		names(ArcList) <- c( "sid", "recid", "bff", "wid" ) # adds names to columns
		ArcList <- with( ArcList, ArcList[ order( sid, recid, wid), ] )

### Create a separate set of records for each wave
### (code needs to be customized for each selected wave set)
### Select by wave, only for recievers where nominee is "best friend" (bff)
### (NOTE: does not account for 0's, i.e. no bff's chosen. Done later.)

		SAff.1 <- with(ArcList, ArcList[ wid == 1, ] ) #extracts edges in wave 1
		SAff.2 <- with(ArcList, ArcList[ wid == 2, ] ) #extracts edges in wave 2
		SAff.3 <- with(ArcList, ArcList[ wid == 3, ] ) #extracts edges in wave 3
		n <- 50 # this is the number of nodes which is not provided by the arclist

### and has to be repeated separately for each of the waves
### (you may loop over the waves if you like),

### For transforming a matrix into an adjacency list

## create indicator matrix of non-zero entries of a
        ones <- !friend.data.w1 %in% 0
## create empty edge list of desired length
        edges <- matrix(0, sum(ones), 3)
# fill the columns of the edge list
        edges[, 1] <- row(friend.data.w1)[ones]
        edges[, 2] <- col(friend.data.w1)[ones]
        edges[, 3] <- friend.data.w1[ones]
# if desired, order edge list by senders and then receivers
        edges <- edges[order(edges[, 1], edges[, 2]), ]

### For transforming an arclist into a matrix
	## First remove the fourth coulmn indicating the wave, so that we are left
	## with sender, receiver and value of the tie, and make it into matrix format
		SAff.1.copy <- SAff.1[, 1:3]
		SAff.1.copy <- as.matrix(SAff.1.copy)
# create empty adjacency matrix
        adj <- matrix(0, 50, 50)
# put edge values in desired places
        adj[edges[, 1:2]] <- edges[, 3]

###################################################################################
################ - READING IN PAJEK DATA - ########################################
###
### If you have data in Pajek format you can use the package "network" in order
### to convert it to a network object. This example is from ?read.paj
###   require( network )
###
###   par( mfrow = c( 2, 2 ) )
###
###test.net.1 <- read.paj("http://vlado.fmf.uni-lj.si/pub/networks/data/GD/gd98/A98.net" )
###   plot( test.net.1,main = test.net.1$gal$title )
###   test.net.2 <- read.paj("http://vlado.fmf.uni-lj.si/pub/networks/data/mix/USAir97.net" )
###   plot( test.net.2,main = test.net.2$gal$title )
###
###################################################################################

# Before we work with the data, we want to be sure it is correct. A simple way
# to check that our data is a matrix is the command class()

        class( friend.data.w1 )

# to list the properties of an object attributes( friend.data.w1 )
# (different classes have different attributes)

# To check that all the data has been read in, we can use the dim() command.
# The adjacency matrix should have the same dimensions as the original data
# (here, 50 by 50).

        dim(friend.data.w1)
        dim(drink)

# To check the values are correct, including missing values, we can use
# the following commands to tabulate the variables.

        table( friend.data.w1, useNA = 'always' )
        table( friend.data.w2, useNA = 'always' )
        table( friend.data.w3, useNA = 'always' )
        table( drink, useNA = 'always' )
        table( smoke, useNA = 'always' )

# NA is the R code for missing data (Not Available).
# This data set happens to have no missings
# (see the data description on the Siena website).
# If there are any missings, it is necessary to tell R about the missing data codes.
# Let us do as if the missing codes for the friendship network were 6 and 9.
# This leads to the following commands.
# (For new R users: the c() function used here as "c(6,9)" constructs
#  a vector [c for concatenate] consisting of the numbers 6 and 9.
#  This function is used a lot in basic R.)

        friend.data.w1[ friend.data.w1 %in% c(6,9) ] <- NA
        friend.data.w1[ friend.data.w2 %in% c(6,9) ] <- NA
        friend.data.w1[ friend.data.w3 %in% c(6,9) ] <- NA

# Commands for descriptive analysis are in the script: RSienaSNADescriptives.R

############## - SELECTING SUBSETS OF DATA - ###################################

# To select a subset of the data based on an actor variable, say,
# those who have the value 2 or 3 on drinking at time 1
# (the possibilities are endless, but hopefully this will serve as a pattern)

      use <- drink[, 1] %in% c(2, 3)

# This creates a logical vector which is TRUE for the cases where the condition
# is satisfied. To view or check, display the vectors:

     	drink[ , 1 ]
     	use

# and the number of selected cases is displayed by

      sum( use )

# or

	table( use )

# To have this arrayed more neatly side by side, you can create and display
# a matrix with the desired information,

        aa <- matrix(nrow=50, ncol=2)
        aa[,1] <- drink[,1]
        aa[,2] <- use
        aa
	#the first column contains the drink level and the second whether this
	#level is 2 or 3 (1) or not (0)

# or as a shorter alternative

	aa <- cbind(drink[ , 1 ],use)
	aa

# Given this selection, submatrices can be formed in case the analyses
# are to be done for this subset only:

        friend1.data.w1 <- friend.data.w1[ use, use ]
        friend1.data.w2 <- friend.data.w2[ use, use ]
        drink1 <- drink[ use, ]



################################################################################
###
### ---- PROCEED TO RscriptSNADescriptives.R FOR DESCRIPTIVE ANAYLSIS ----------
###
################################################################################
