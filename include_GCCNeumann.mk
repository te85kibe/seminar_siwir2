CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

CFLAGS   = -O3 -Wall -Winline -Wshadow -std=c++11 -DNEUMANN
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
CPPFLAGS = $(CFLAGS)
LFLAGS   =  
DEFINES  =
INCLUDES =
LIBS     =
