CC  = gcc
CXX = g++
FC  = gfortran
LINKER = $(CXX)

CFLAGS   = -O3 -Wno-format  -Wall -DNDEBUG -Winline -Wshadow -std=c++11
CXXFLAGS = $(CFLAGS)
FCFLAGS  = 
CPPFLAGS = -std=c++0x
LFLAGS   =  
DEFINES  = -D_GNU_SOURCE
INCLUDES =
LIBS     =
