# This file is part of t8dg
# Non-recursive Makefile.am in exec
# Included from toplevel directory

bin_PROGRAMS += \
	exec/t8dg 

exec_t8dg_SOURCES = \
	thirdparty/mptrac/libtrac.c \
	src/main.c

exec_t8dg_LDADD = $(LDADD) thirdparty/mptrac/libt8dgmptrac.la

