
# This file is part of t8dg
# Makefile.am in toplevel directory

# include makefiles from installed subpackages
ACLOCAL_AMFLAGS = -I config @T8DG_SC_AMFLAGS@ @T8DG_P4EST_AMFLAGS@ @T8DG_T8_AMFLAGS@
if T8DG_SC_MK_USE
@T8DG_SC_MK_INCLUDE@
endif
if T8DG_P4EST_MK_USE
@T8DG_P4EST_MK_INCLUDE@
endif
if T8DG_T8_MK_USE
@T8DG_T8_MK_INCLUDE@
endif

# initialize empty variables
AM_CPPFLAGS =
CLEANFILES =
DISTCLEANFILES =
EXTRA_DIST =
LDADD =
LINT_CSOURCES =
TESTS =
bin_PROGRAMS =
check_PROGRAMS =
include_HEADERS =
lib_LTLIBRARIES =
nodist_include_HEADERS =
nobase_dist_include_HEADERS =
noinst_HEADERS =
noinst_PROGRAMS =
sysconf_DATA =
dist_t8dgdata_DATA =

# use this if you want to link in T8DG without autotools
sysconf_DATA += Makefile.t8dg.mk
CLEANFILES += Makefile.t8dg.mk
Makefile.t8dg.mk : Makefile.t8dg.pre
	cat $< | \
        sed -e 's,{\(\(.*prefix\|sysconfdir\)\)},{t8dg_\1},g' \
            -e 's,^\(\(.*prefix\|sysconfdir\) *=\),t8dg_\1,g' > $@

# install T8DG m4 macros in the correct directory
t8dgaclocaldir = $(datadir)/aclocal
dist_t8dgaclocal_DATA = config/t8dg_include.m4

# install t8dg data in the correct directory
t8dgdatadir = $(datadir)/data

# setup test environment
if T8DG_MPIRUN
LOG_COMPILER = @T8DG_MPIRUN@
AM_LOG_FLAGS = @T8DG_MPI_TEST_FLAGS@
endif

# recurse only into subpackages
SUBDIRS = @T8DG_SC_SUBDIR@ @T8DG_P4EST_SUBDIR@ @T8DG_T8_SUBDIR@
DIST_SUBDIRS = $(SUBDIRS)

# handle toplevel directory
EXTRA_DIST += \
  bootstrap build-aux/git-version-gen build-aux/git2cl doc
DISTCLEANFILES += \
        _configs.sed src/t8dg_config.h @T8DG_DISTCLEAN@
.PHONY: ChangeLog

## # Force libtool to link with c++
## nodist_EXTRA_T8DG_la_SOURCES = dummy.cxx
## sourcescheck_DUMMY = dummy.cxx

# non-recursive build
include src/Makefile.am
include exec/Makefile.am
include test/Makefile.am
include thirdparty/mptrac/Makefile.am

# revision control and ChangeLog
ChangeLog:
	(GIT_DIR=@top_srcdir@/.git git log > .ChangeLog.tmp && \
         cat .ChangeLog.tmp | @top_srcdir@/build-aux/git2cl > ChangeLog) ; \
        rm -f .ChangeLog.tmp

dist-hook:
	echo $(VERSION) > $(distdir)/.tarball-version
	test "x$(VERSION)" = "x`@top_srcdir@/build-aux/git-version-gen \
          @top_srcdir@/.tarball-version`" || \
        ((echo "Stale version;" ; echo "Please run:" ; \
          echo "     cd @top_srcdir@ && ./bootstrap" ; \
          echo "before make dist") 1>&2 ; rm -r $(distdir) ; exit 1)
if T8DG_DIST_DENY
	@echo "-----------------------------------------------------"
	@echo "make dist does not work with external subpackages"
	@echo "-----------------------------------------------------"
	rm -rf $(distdir)
	@exit 1
endif

clean-local:
	rm -f ChangeLog *vtu *.visit

maintainer-clean-local:
	rm -rf doxygen

doxygen: Doxyfile
	doxygen

.PHONY: doxygen
