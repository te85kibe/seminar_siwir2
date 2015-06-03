CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

CFLAGS   = -O3 -Wall -Winline -Wshadow -fopenmp -std=c++11 
CXXFLAGS = $(CFLAGS) 
FCFLAGS  = $(CFLAGS) 
CPPFLAGS = $(CFLAGS)
LFLAGS   = $(CFLAGS) 
DEFINES  = 
INCLUDES =
LIBS     =
