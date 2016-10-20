#/******************************************************************
# * Modul name : makefile
# * Author     : Zoltan Kato
# * Description:
# * Make file for making ColorMRFdemo
# * make all             = compilation + installation + clean.
# * make compile         = compilation.
# * make install         = compilation + installation.
# * make clean           = clean.
# * 
#
# Adapted from the wxwindows sample makefile written by Robert Roebling. 
#
# This makefile requires a Unix version of wxWidgets
# to be installed on your system. 
#
#
# * $Id: makefile,v 1.1 2009/01/09 19:59:12 kato Exp $
# * $Revision: 1.1 $
# * $State: Exp $
# * $Log: makefile,v $
# * Revision 1.1  2009/01/09 19:59:12  kato
# * Initial revision
# *
# *
# *****************************************************************/
#
# dirs, flags & libs
#

OS=$(shell uname -s)
CFLAGS=$(shell wx-config --cflags) -D$(OS) -D__USE_FIXED_PROTOTYPES__
CFLAGS+= -Wno-potentially-evaluated-expression -Wno-comment

LDFLAGS= $(shell wx-config --libs)
LIBS= -lm
ifeq ($(OS),SunOS)
LIBS=$(LIBS) -liberty
endif
VPATH= ../src

#
# Source files:
#
CFILES= $(wildcard ../src/*.cpp)
#
# Target files
#
OBJECTS=$(subst ../src/,./,$(CFILES:.cpp=.o))
TARGETS= colormrfdemo
#
# Suffixes
#
.SUFFIXES: .cpp.o
.cpp.o:
	g++  $(CFLAGS) -c $^
#
# Rules
#

$(TARGETS) : $(OBJECTS)
	g++ -o $@ $^ $(LDFLAGS) $(LIBS)

compile: $(TARGETS)

clean:
	rm -f $(OBJECTS)
	rm -f $(TARGETS)
