# This file is part of t8dg
# Non-recursive Makefile.am in exec
# Included from toplevel directory

bin_PROGRAMS += \
	exec/t8dg \
	exec/t8dg_try_with_petsc

exec_t8dg_SOURCES = \
	src/main.c

exec_t8dg_try_with_petsc_SOURCES = \
	exec/t8dg_try_with_petsc.c

exec_t8dg_try_with_petsc_CPPFLAGS = @T8DG_PETSC_INCLUDE_DIRS@ $(AM_CPPFLAGS)
exec_t8dg_try_with_petsc_LDADD = @T8DG_PETSC_LINK_LIBS@ $(LDADD)