# This file is part of t8dg
# Non-recursive Makefile.am in test
# Included from toplevel directory

t8dg_googletest_programs = \
  test/t8dg_tests

test_t8dg_tests_SOURCES = test/t8dg_test_main.cxx \
  test/t8dg_test_dmatrix.cxx \
  test/t8dg_test_functionbasis.cxx \
  test/t8dg_test_quadrature.cxx \
  test/t8dg_test_timestepping.cxx \
  test/t8dg_test_values.cxx
test_t8dg_tests_LDADD = $(LDADD) test/libgtest.la
test_t8dg_tests_LDFLAGS = $(AM_LDFLAGS) -pthread
test_t8dg_tests_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/thirdparty/googletest/ \
        -pthread \
  -std=c++11

# Build Googletest library
check_LTLIBRARIES = test/libgtest.la
test_libgtest_la_SOURCES = thirdparty/googletest/gtest/gtest-all.cc
test_libgtest_la_CPPFLAGS = -I$(top_srcdir)/thirdparty/googletest -std=c++11
test_libgtest_la_LDFLAGS = -pthread


TESTS +=  $(t8dg_googletest_programs)
check_PROGRAMS += $(t8dg_googletest_programs)
