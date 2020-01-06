# RSiena new version 1.2-21 (2019-12-18)

## The only remaining NOTE was that on x86_64-pc-linux-gnu (64-bit)
   apparently the file Siena.out was not deleted. This now is done
   by calling unlink().

# RSiena new version 1.2-20 (2019-12-16)

## The note in Debian about existence of file Siena.out in the 
   check directory was fixed, using file.remove().

## The additional issue for the earlier version 1.2-12 in LTO 
   (about getChainProbabilities) was fixed by dropping this function 
   which was not used anyway.

## The additional issue for the earlier version 1.2-12 in rchk (about DumpChain)
   was fixed by dropping this function which was not used anyway.


# RSiena new version 1.2-19 (2019-12-16)

## Kurt Hornik sent me a message that it was necessary to replace the use of
   (class(..) == ... ) by (inherits(...)). 
   This was done, and fixes the error in the earlier version 1.2-12 occurring 
   for r-devel-linux-x86_64-debian-clang.

## Passed checks on Windows, Mac, and Linux for R-release and R-devel.
* No ERRORs or WARNINGs.
* On some systems there is a NOTE about the installed package size,
  which is due to the use of a lot of compiled C++ code.
* Sometimes NOTEs about computation times of examples.
  These require simulations, and times are barely over 5s. 

## There are quite some parts in the .Rd files with donttest.
   This is because executing them is too time-consuming.
* R CMD check --run-donttest on Windows OK.

## Reverse dependencies checked: btergm.

## Reverse suggests checked: netdiffuse.