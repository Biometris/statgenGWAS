## require at least autoconf 2.61
AC_PREREQ(2.61)

## Process this file with autoconf to produce a configure script.
AC_INIT(statgenGWAS, m4_esyscmd_s([awk '/^Version:/ {print $2}' DESCRIPTION]))

## Set R_HOME, respecting an environment variable if one is set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi
## Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

AC_OUTPUT
