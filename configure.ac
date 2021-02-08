#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ
AC_INIT(bin/LDBlockShow, 1.36, hewm2008@gmail.com)
AM_INIT_AUTOMAKE(bin/LDBlockShow,1.36)
AC_CONFIG_SRCDIR([src/LDBlockShow.cpp])
AC_CONFIG_HEADER([config.h])
# Checks for programs
AC_PROG_AWK
AC_PROG_CXX([g++])
AC_PROG_CC([gcc])
AC_PROG_CPP
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lz':
AC_CHECK_LIB([z], [main],,AC_MSG_ERROR([You need zlib >= 1.2.3 to build $PACKAGE]))


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([zlib.h],[],[AC_MSG_ERROR([You need include zlib.h to build $PACKAGE])],[])
AC_CHECK_HEADERS([math.h],[],[AC_MSG_ERROR([You need include math.h to build $PACKAGE])],[])


# Checks for typedefs, structures, and compiler characteristics.

#AC_CHECK_HEADER_STDBOOL
#AC_C_INLINE
#AC_TYPE_INT32_T
#AC_TYPE_SIZE_T
#AC_TYPE_UINT32_T
#AC_TYPE_UINT8_T



# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset regcomp sqrt strdup strncasecmp strtol])




AC_CONFIG_FILES([Makefile])
AC_OUTPUT
