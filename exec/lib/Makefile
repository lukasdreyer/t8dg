# Googletest source folder
GTEST_DIR=../../thirdparty/googletest

# C++ Compiler and flags
CXX=mpicxx
CXXFLAGS=-c --std=c++11 -O3
INCLUDE=-I${GTEST_DIR}

.PHONY: clean all

all: googletest

gtest-all.o:
	${CXX} ${CXXFLAGS} ${INCLUDE} ${GTEST_DIR}/gtest/gtest-all.cc -o gtest-all.o

libgtest.a: gtest-all.o
	ar rcs libgtest.a gtest-all.o

googletest: libgtest.a

clean:
	rm gtest-all.o libgtest.a
