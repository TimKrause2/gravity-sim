CPPFLAGS=-g
LDLIBS=-lm
CC=g++

gsim:gsim.o


gsim.o:gsim.cpp vector.h gsim.h
