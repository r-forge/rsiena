# RSiena new version 1.2-19 (2019-12-16)

## Kurt Hornik sent me a message that it was necessary to replace the use of
   (class(..) == ... ) by (inherits(...)). This was done.

## Passed checks on Windows, Mac, and Linux for R-release and R-devel.
* No ERRORs or WARNINGs.
* On some systems there is a NOTE about the installed package size,
  which is due to the use of a lot of compiled C++ code.

## There are quite some parts in the .Rd files with donttest.
   This is because executing them is too time-consuming.
* R CMD check --run-donttest on Windows OK.

## Reverse dependencies checked: btergm.

## Reverse suggests checked: netdiffuse.