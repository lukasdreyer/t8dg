dnl
dnl This file is part of t8dg.
dnl

AC_INIT([t8dg],
        [m4_esyscmd([build-aux/git-version-gen .tarball-version])],
        [johannes.holke@dlr.de])
AC_PREREQ(2.61)
AC_CONFIG_HEADERS([src/pre_config.h])
AC_CONFIG_SRCDIR([src/t8dg.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])
AX_PREFIX_CONFIG_H([src/t8dg_config.h])
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])
AM_SILENT_RULES
SC_VERSION([T8DG])

echo "o---------------------------------------"
echo "| Checking options"
echo "o---------------------------------------"

T8DG_ARG_ENABLE([debug], [enable debug mode (assertions and extra checks)],
                 [DEBUG])

echo "o---------------------------------------"
echo "| Checking MPI and related programs"
echo "o---------------------------------------"

dnl A nonempty second/third argument causes to enable F77+F90/CXX, respectively.
SC_MPI_CONFIG([T8DG], [yes], [yes])
SC_MPI_ENGAGE([T8DG])
# This is needed for compatibility with automake >= 1.12
m4_ifdef([AM_PROG_AR],[AM_PROG_AR])
SC_PROG_LINT
SC_C_VERSION
LT_INIT

echo "o---------------------------------------"
echo "| Checking libraries"
echo "o---------------------------------------"

SC_CHECK_LIBRARIES([T8DG])
P4EST_CHECK_LIBRARIES([T8DG])
T8_CHECK_LIBRARIES([T8DG])
T8DG_CHECK_LIBRARIES([T8DG])

echo "o---------------------------------------"
echo "| Checking headers"
echo "o---------------------------------------"

dnl AC_CHECK_HEADERS([arpa/inet.h netinet/in.h unistd.h])

echo "o---------------------------------------"
echo "| Checking functions"
echo "o---------------------------------------"

dnl AC_CHECK_FUNCS([fsync])

echo "o---------------------------------------"
echo "| Checking subpackages"
echo "o---------------------------------------"

T8DG_DIST_DENY=
T8DG_DISTCLEAN=
SC_AS_SUBPACKAGE([T8DG])
export P4EST_SC_SOURCE=../sc
export T8_SC_SOURCE=../sc
export T8_P4EST_SOURCE=../p4est
P4EST_AS_SUBPACKAGE([T8DG])
T8_AS_SUBPACKAGE([T8DG])
AM_CONDITIONAL([T8DG_DIST_DENY], [test "x$T8DG_DIST_DENY" != x])
AC_SUBST([T8DG_DISTCLEAN])

# Print summary.

AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C preprocessor flags])
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C compiler flags])
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

echo "o----------------------------------"
echo "| Results for t8dg are"
echo "o----------------------------------"
echo "| CPP:         $CPP"
echo "| CPPFLAGS:    $CPPFLAGS"
echo "| CC:          $CC"
echo "| CFLAGS:      $CFLAGS"
echo "| LDFLAGS:     $LDFLAGS"
echo "| LIBS:        $LIBS"
echo "o----------------------------------"

# Create output files.
AC_CONFIG_FILES([Makefile Makefile.t8dg.pre Doxyfile])
AC_OUTPUT

# Final messages.
SC_FINAL_MESSAGES([T8DG])
P4EST_FINAL_MESSAGES([T8DG])
T8_FINAL_MESSAGES([T8DG])
T8DG_FINAL_MESSAGES([T8DG])
