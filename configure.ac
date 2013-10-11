#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

m4_define([XXVERSION], m4_esyscmd([grep VERSION_STRING version.h | sed 's/^[^"]*"\([^"\\]*\).*/\1/' | tr -d '\n']))

AC_INIT([libsoxconvolver], [XXVERSION], [honeycomb77@gmail.com])
AC_CONFIG_SRCDIR([dft_filter.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_PROG_INSTALL
AM_PROG_CC_C_O
LT_LANG([Windows Resource])
AM_CONDITIONAL([HAVE_WINDRES],
               [test x$ac_cv_prog_ac_ct_RC != x -o x$ac_cv_prog_RC != x])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([memalign _aligned_malloc])

VI_CURRENT=0
VI_AGE=0
VI_REVISION=0
AC_SUBST(VI_CURRENT)
AC_SUBST(VI_AGE)
AC_SUBST(VI_REVISION)

AC_CANONICAL_HOST

case ${host_cpu} in
i?86) XX_IX86="yes";;
esac

case ${CC} in
*gcc*) XX_GCC="yes";;
*g++*) XX_GCC="yes";;
esac

AM_CONDITIONAL([IX86],[test "x${XX_IX86}" != "x" ])
AM_CONDITIONAL([GCC],[test "x${XX_GCC}" != "x" ])

AC_ARG_ENABLE([simd],
              [AS_HELP_STRING([--disable-simd],
                              [Disable SIMD (default is enabled)])],
              [],
              [enable_simd=yes])
AS_IF([test "x$enable_simd" == xno],
      [AC_DEFINE([DISABLE_SIMD], [1], [Disable SIMD support])])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
