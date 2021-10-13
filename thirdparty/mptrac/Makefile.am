# This file is part of t8dgcode
# Non-recursive Makefile.am in example/mptrac/thirdparty/mptrac
# Included from example/mptrac directory

libt8dgmptrac_installed_headers = \
  thirdparty/mptrac/libtrac.h

libt8dgmptrac_compiled_sources = \
  thirdparty/mptrac/libtrac.c

lib_LTLIBRARIES += thirdparty/mptrac/libt8dgmptrac.la
thirdparty_mptrac_libt8dgmptrac_la_SOURCES = \
	$(libt8dgmptrac_internal_headers) \
	$(libt8dgmptrac_compiled_sources)
thirdparty_mptrac_libt8dgmptrac_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8DG_CPPFLAGS)
thirdparty_mptrac_libt8dgmptrac_la_LDFLAGS =
EXTRA_thirdparty_mptrac_libt8dgmptrac_la_DEPENDENCIES =

nodist_include_HEADERS += $(libt8dgmptrac_generated_headers)
nobase_dist_include_HEADERS += $(libt8dgmptrac_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/../thirdparty/mptrac
