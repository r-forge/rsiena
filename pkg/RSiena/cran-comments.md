# RSiena new version 1.2-3 (2017-9-08)

## Issue for previous version (1.1-232) on gcc-UBSAN (runtime error) solved.

## Warning for previous version (1.1-232) on windows-devel
(cleanup: Non-Windows OSes require LF line endings) solved.

##Test environments
* local Windows 7, R 3.4.1: OK
* local Windows 7, R-devel (2017-8-24) through devtools: OK
* local Windows 7, R-patched (2017-8-25) through devtools: OK
* R-Hub Debian Linux, R-devel, GCC: OK
* R-Hub macOS 10.9 Mavericks, R-oldrel (experimental): OK
* R-Hub macOS 10.11 El Capitan, R-release (experimental): OK
* R-Hub Ubuntu Linux 16.04 LTS, R-devel, GCC: OK

## R CMD check results
* All results were OK, no ERRORs, WARNINGs, or NOTEs.