# This file is part of t8dgcode
# Non-recursive Makefile.am in src
# Included from toplevel directory

libt8dg_generated_headers = src/t8dg_config.h
libt8dg_installed_headers = \
  src/t8dg.h \
  src/t8dg_functionbasis.h \
  src/t8dg_refcount.h \
  src/t8dg_global_values.h \
  src/t8dg_local_values.h \
  src/t8dg_coarse_geometry.h \
  src/t8dg_values.h \
  src/t8dg_vertexset.h \
  src/t8dg_flux.h \
  src/t8dg_flux_implementation.h \
  src/t8dg_flux_data_base.h \
  src/t8dg_dof.h \
  src/t8dg_quad.h \
  src/t8dg_mortar.h \
  src/t8dg_mptrac.h \
  src/t8dg_timestepping.h \
  src/t8dg_dmatrix.h \
  src/t8dg_quadrature.h \
  src/t8dg_sc_array.h \
  src/t8dg_tensor.h \
  src/t8dg_geometry.h \
  src/t8dg_common.h \
  src/t8dg_adapt.h \
  src/t8dg_output.h \
  src/t8dg_cmesh.h \
  src/t8dg_advect_diff_problem.h

libt8dg_internal_headers =

libt8dg_compiled_sources = \
  src/t8dg_functionbasis.c \
  src/t8dg_refcount.c \
  src/t8dg_global_values.c \
  src/t8dg_sc_array.c \
  src/t8dg.cxx \
  src/t8dg_local_values.cxx \
  src/t8dg_tensor.c \
  src/t8dg_coarse_geometry.c \
  src/t8dg_mortar.cxx \
  src/t8dg_mptrac.cxx \
  src/t8dg_timestepping.c \
  src/t8dg_dmatrix.c \
  src/t8dg_values.cxx \
  src/t8dg_dof.c \
  src/t8dg_quad.c \
  src/t8dg_vertexset.c \
  src/t8dg_flux.cxx \
  src/t8dg_flux_implementation.cxx \
  src/t8dg_quadrature.c \
  src/t8dg_geometry.cxx \
  src/t8dg_common.c \
  src/t8dg_adapt.cxx \
  src/t8dg_output.cxx \
  src/t8dg_cmesh.cxx \
  src/t8dg_advect_diff_problem.cxx

# this variable is used for headers that are not publicly installed
T8DG_CPPFLAGS =

lib_LTLIBRARIES += src/libt8dg.la
src_libt8dg_la_SOURCES = \
	$(libt8dg_internal_headers) \
	$(libt8dg_compiled_sources)
src_libt8dg_la_CPPFLAGS = $(AM_CPPFLAGS) $(T8DG_CPPFLAGS)
src_libt8dg_la_LDFLAGS = -release $(VERSION)
LDADD += @top_builddir@/src/libt8dg.la \
        @T8DG_T8_LDADD@ @T8DG_P4EST_LDADD@ @T8DG_SC_LDADD@
EXTRA_src_libt8dg_la_DEPENDENCIES = @T8DG_T8_LDADD@ @T8DG_P4EST_LDADD@ @T8DG_SC_LDADD@

nodist_include_HEADERS += $(libt8dg_generated_headers)
include_HEADERS += $(libt8dg_installed_headers)

AM_CPPFLAGS += -I@top_srcdir@/src \
        @T8DG_SC_CPPFLAGS@ @T8DG_P4EST_CPPFLAGS@ @T8DG_T8_CPPFLAGS@
