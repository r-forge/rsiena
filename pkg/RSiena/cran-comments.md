# RSiena new version 1.2-23 (2020-01-12)

## Passed checks on Windows, Mac, and Linux for R-release and R-devel.
* No ERRORs or WARNINGs.
* On some systems there is a NOTE about the installed package size,
  which is due to the use of a lot of compiled C++ code.
* Sometimes NOTEs about computation times of examples.
  These require simulations, and times are barely over 5s. 

## Kurt Hornik sent me a message that it was necessary to replace the use of
   (class(..) == ... ) by (inherits(...)). 
   This was done, and fixes the error in the earlier version 1.2-12 occurring 
   for r-devel-linux-x86_64-debian-clang.

## Brian Ripley sent me a message that there may be minor changes in output
   depending on the platform. This now is avoided by minor changes in the tests.
   Brian also sent me a message about strange things in the Configure files.
   Since tests are passed (as above) on various R-hub platforms
   I think this is not a major obstacle for now, and will look into it later.

## There are quite some parts in the .Rd files with donttest.
   This is because executing them is too time-consuming.
* R CMD check --run-donttest on Windows OK.

## Reverse dependencies checked: btergm.

## Reverse suggests checked: netdiffuse.