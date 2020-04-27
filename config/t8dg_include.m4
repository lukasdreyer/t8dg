dnl
dnl t8dg_include.m4 - custom macros
dnl

dnl Documentation for macro names: brackets indicate optional arguments

dnl T8DG_ARG_ENABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable T8DG_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional T8DG_TOKEN
dnl Default is disabled
dnl
AC_DEFUN([T8DG_ARG_ENABLE],
         [SC_ARG_ENABLE_PREFIX([$1], [$2], [$3], [T8DG])])

dnl T8DG_ARG_DISABLE(NAME, COMMENT, TOKEN)
dnl Check for --enable/disable-NAME using shell variable T8DG_ENABLE_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If enabled, define TOKEN to 1 and set conditional T8DG_TOKEN
dnl Default is enabled
dnl
AC_DEFUN([T8DG_ARG_DISABLE],
         [SC_ARG_DISABLE_PREFIX([$1], [$2], [$3], [T8DG])])

dnl T8DG_ARG_WITH(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable T8DG_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional T8DG_TOKEN
dnl Default is without
dnl
AC_DEFUN([T8DG_ARG_WITH],
         [SC_ARG_WITH_PREFIX([$1], [$2], [$3], [T8DG])])

dnl T8DG_ARG_WITHOUT(NAME, COMMENT, TOKEN)
dnl Check for --with/without-NAME using shell variable T8DG_WITH_TOKEN
dnl If shell variable is set beforehand it overrides the option
dnl If with, define TOKEN to 1 and set conditional T8DG_TOKEN
dnl Default is with
dnl
AC_DEFUN([T8DG_ARG_WITHOUT],
         [SC_ARG_WITHOUT_PREFIX([$1], [$2], [$3], [T8DG])])

dnl T8DG_CHECK_LIBRARIES(PREFIX)
dnl This macro bundles the checks for all libraries and link tests
dnl that are required by T8DG.  It can be used by other packages that
dnl link to T8DG to add appropriate options to LIBS.
dnl
AC_DEFUN([T8DG_CHECK_LIBRARIES],
[
T8DG_CHECK_CPPSTD([$1])
])

dnl T8DG_AS_SUBPACKAGE(PREFIX)
dnl Call from a package that is using T8DG as a subpackage.
dnl Sets PREFIX_DIST_DENY=yes if T8DG is make install'd.
dnl
AC_DEFUN([T8DG_AS_SUBPACKAGE],
         [SC_ME_AS_SUBPACKAGE([$1], [m4_tolower([$1])], [T8DG], [t8dg])])

dnl T8DG_FINAL_MESSAGES(PREFIX)
dnl This macro prints messages at the end of the configure run.
dnl
AC_DEFUN([T8DG_FINAL_MESSAGES],
[
])
