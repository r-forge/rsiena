# Build System: See doc/RSienaDeveloper.tex for explanations.

# High level targets for this file:
#
# clean
#   Remove temporary files.
# build, check, install
#   Correspond the to `R CMD *` versions.
# test_capture, test
#   Trigger the regression testing. regsave captures the state, regtest
#   tests against the saved state.

# Parse DESCRIPTION file.
PKG_VERSION := $(shell grep -i ^version DESCRIPTION | cut -d\  -f2-)
PKG_REVISION := $(shell echo $(PKG_VERSION) | cut -d- -f2-)
PKG_NAME := $(shell grep -i ^package DESCRIPTION | cut -d\  -f2-)
PKG_DATE := $(shell grep -i ^date DESCRIPTION | cut -d\  -f2-)
PKG_TARBALL := $(PKG_NAME)_$(PKG_VERSION).tar.gz
PKG_IMPORTS := $(shell grep -i ^imports DESCRIPTION | cut -d\  -f2-)
PKG_SUGGESTS := $(shell grep -i ^suggests DESCRIPTION | cut -d\  -f2-)

# Actual source information.
# SRC_REVISION := $(shell svn info --show-item revision)
SRC_LIST = src/sources.list

# Top level source files containing the R binding.
SRC_BINDING = \
	RInterface.cpp \
	RUtil.cpp \
	siena07internals.cpp \
	siena07models.cpp \
	siena07setup.cpp \
	siena07utilities.cpp

# Folders in src/ needed entirely.
SRC_MODULES = data Eigen estimator logger model network sim utils

# R command settings.
R = R
R_RUN = $(R) -e
R_BUILD = $(R) CMD build
R_CHECK = $(R) CMD check --as-cran
R_INSTALL = $(R) CMD INSTALL

.PHONY:
all: clean check

.PHONY:
clean:
# ifeq ($(OS),Windows_NT)
# 	sh ./cleanup.win
# else
# 	./cleanup
# endif
	-rm -rf $(PKG_NAME).Rcheck
	-rm -f $(PKG_TARBALL)

.PHONY:
check: $(PKG_TARBALL)
	$(R_CHECK) $<

.PHONY:
install: $(PKG_TARBALL)
	$(R_INSTALL) $<

.PHONY:
build: $(PKG_TARBALL)
$(PKG_TARBALL): configure $(SRC_LIST)
	$(R_BUILD) .

.PHONY:
commitchecks: # test
	# test $(PKG_REVISION) -eq $$(($(SRC_REVISION)+1)) # DESCRIPTION revision
	test $(PKG_DATE) = $$(date -I) # DESCRIPTION date
	chmod a+x configure cleanup # script permissions

configure: configure.ac
	autoconf

# build the list of source files
$(SRC_LIST): $(shell find src -iname '*.cpp')
	echo -n "$(SRC_BINDING)" >$@
	cd src && find $(SRC_MODULES) -iname '*.cpp' -printf ' %p' >>$(@:src/%=%)

# dependencies
.PHONY:
install_dep:
	$(R_RUN) "options(repos='http://stat.ethz.ch/CRAN');\
	pkgs <- c(strsplit('$(PKG_IMPORTS)', ', '), strsplit('$(PKG_SUGGESTS)', ', '));\
	for (pkg in pkgs) install.packages(pkg)"

# regression tests
.PHONY:
test_capture:
	cd inst/unitTests && \
	$(R_RUN) "require($(PKG_NAME)); record_values <- T; RSienaTest:::run_tests(dir='.')"

.PHONY:
test:
	cd inst/unitTests && \
	$(R_RUN) "require($(PKG_NAME)); RSienaTest:::run_tests(dir='.')"

.PHONY:
memtest:
	R CMD INSTALL .
	R -d 'valgrind --leak-check=full' -f tests/parallel.R
