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
INCLUDE=$(HOME)/include
LIBPATH=$(HOME)/lib-$(ARCH)
INSTALLDIR=$(HOME)/bin-$(ARCH)
MANDIR=$(HOME)/man
OS=`uname -s`
CFLAGS= `wx-config --cflags` -D$(OS) -I$(INCLUDE) -ansi -D__USE_FIXED_PROTOTYPES__
LDFLAGS= -L$(LIBPATH) `wx-config --libs`
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
.cpp.o:	; echo 'Compiling $*.cpp'; gcc  $(CFLAGS) -c $^
#
# Rules
#
.SILENT:

all: install clean

$(TARGETS) : $(OBJECTS)
	echo 'Building $@'
	gcc -o $@ $^ $(LDFLAGS) $(LIBS)

compile: $(TARGETS)

clean:
	echo 'Cleaning up'
	/bin/rm -f $(OBJECTS)
	/bin/rm -f $(TARGETS)

install: compile
	echo 'Installing $(TARGETS) in $(INSTALLDIR)/'
	cp $(TARGETS) $(INSTALLDIR)
