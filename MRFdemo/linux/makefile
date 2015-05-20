#/******************************************************************
# * Modul name : makefile
# * Author     : Zoltan Kato
# * Description:
# * Make file for making MRFdemo
# * make all             = compilation + installation + clean.
# * make compile         = compilation.
# * make install         = compilation + installation.
# * make clean           = clean.
# * 
#
# Adapted from the wxwindows sample makefile written by Robert Roebling. 
#
# This makefile requires a Unix version of wxWindows
# to be installed on your system. This is most often
# done typing "make install" when using the complete
# sources of wxWindows or by installing the two
# RPM packages wxGTK.XXX.rpm and wxGTK-devel.XXX.rpm
# under Linux.
#
#
# * $Id: makefile,v 1.2 2005/01/13 11:15:53 kato Exp $
# * $Revision: 1.2 $
# * $State: Exp $
# * $Log: makefile,v $
# * Revision 1.2  2005/01/13 11:15:53  kato
# * updated for new directory structure
# *
# * Revision 1.1  2004/12/08 13:15:24  kato
# * Initial revision
# *
# *****************************************************************/
#
# dirs, flags & libs
#
INCLUDE=$(HOME)/include
LIBPATH=$(HOME)/lib-$(ARCH)
INSTALLDIR=$(HOME)/bin-$(ARCH)
MANDIR=$(HOME)/man
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
TARGETS= mrfdemo
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
