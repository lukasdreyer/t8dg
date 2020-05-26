dnl T8DG_CHECK_CBLAS
dnl Check for cblas support and link a test program
dnl
dnl This macro tries to link to the cblas library.
dnl Use the LIBS variable on the configure line to specify a different library
dnl or use --with-cblas=<LIBRARY>
dnl
dnl Using --with-cblas without any argument defaults to -lcblas.
dnl
AC_DEFUN([T8DG_CHECK_CBLAS], [

dnl This link test changes the LIBS variable in place for posterity
dnl SAVE_LIBS="$LIBS"
dnl T8DG_CHECK_LIB([cblas], [nc_open], [CBLAS])
dnl LIBS="$SAVE_LIBS"
dnl AC_MSG_CHECKING([for cblas linkage])

T8DG_ARG_DISABLE([cblas],
  [cblas library (optionally use --enable-cblas=<CBLAS_LIBS>)],
  [CBLAS])
if test "x$T8DG_ENABLE_CBLAS" != xno ; then
  T8DG_CBLAS_LIBS="-lblas"
  if test "x$T8DG_ENABLE_CBLAS" != xyes ; then
    T8DG_CBLAS_LIBS="$T8DG_ENABLE_CBLAS"
    dnl AC_MSG_ERROR([Please provide --enable-cblas without arguments])
  fi
  PRE_CBLAS_LIBS="$LIBS"
  LIBS="$LIBS $T8DG_CBLAS_LIBS"
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
[[
]],[[
  #include <cblas.h>

  cblas_snrm2(0, NULL, 0);
]])],,
                 [AC_MSG_ERROR([Unable to link with cblas library])])
dnl Keep the variables changed as done above
dnl LIBS="$PRE_CBLAS_LIBS"



  AC_MSG_RESULT([successful])
else
  AC_MSG_RESULT([not used])
fi

])

