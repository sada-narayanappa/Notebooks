# Modified for tcp migration project

SHELL = /bin/sh
OS=$(shell uname -s)
FLAVOR=${OS}

ifeq ($(OS),Linux)
	FLAVOR=$(shell cat /etc/os-release | grep ^ID= | cut -d'=' -f2)
else
	ifeq ($(OS),Darwin)
		FLAVOR=${OS}
	endif
endif

CC = g++
PROGRAMS1 = SCORE.EXE_$(FLAVOR)
DEFS = -DHAVE_CONFIG_H -I.
INCLUDES =  -I".."
CPPFLAGS = 
LDFLAGS = 
#LIBS = -lpcap  -lpthread
OBJECTS1 =  CSV.o score.o ncsv.o
DEPENDENCIES = score.o 
CFLAGS = -g -O2 -Wall -Wno-unused
COMPILE = $(CC) -static $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
ifeq ($(FLAVOR),alpine)
	AM_CFLAGS=-static-libstdc++ -static-libgcc
endif
LINK = $(CCLD)  $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

all: $(PROGRAMS1)
.SUFFIXES:
.SUFFIXES: .S .c .o .s

.c.o:
	$(COMPILE) -c $<

.cpp.o:
	$(COMPILE) -c $<

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

$(PROGRAMS1): $(OBJECTS1) $(DEPENDENCIES)
	@rm -f $(PROGRAMS1)
	$(LINK) $(LDFLAGS) $(OBJECTS1) $(LDADD) $(LIBS)


CSV.o: CSV.c CSV.h marray.h
ncsv.o: ncsv.c ncsv.h any.h marray.h
score.o: score.c any.h marray.h


clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS1) SCORE.EXE* core *.core OUT*
