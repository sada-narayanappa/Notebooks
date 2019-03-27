# Modified for tcp migration project
SHELL = /bin/sh

OS=$(shell uname -s)
CC = g++

PROGRAMS = INVX.EXE_$(OS)

DEFS = -DHAVE_CONFIG_H -I.
INCLUDES =  -I".." `pkg-config --cflags eigen3`
CPPFLAGS = 
LDFLAGS = 
LIBS = -lpcap  -lpthread
OBJECTS  =  CSV.o LR.o CSVUtils.o main.o invx.o
DEPENDENCIES = main.o
CFLAGS = -g -O2 -Wall -Wno-unused
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
CLINKS = echo ln -s $@ ${@%.*}

all: $(PROGRAMS) CREATE_LINKS
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

$(PROGRAMS): $(OBJECTS) CSVUtils.o $(DEPENDENCIES)
	@rm -f $(PROGRAMS)
	$(LINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)

$(PROGRAMS1): $(OBJECTS1) $(DEPENDENCIES)
	@rm -f $(PROGRAMS1)
	$(LINK) $(LDFLAGS) $(OBJECTS1) $(LDADD) $(LIBS)

CREATE_LINKS: $(PROGRAMS)
	ln -s  $(PROGRAMS) INVX.EXE

main.o: main.c
CSV.o: CSV.c CSV.h
LR.o: LR.c LR.h
invx.o: invx.c invx.h any.h
CSVUtils.o: CSVUtils.c CSV.h


clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS1) INVX.EXE SCORE.EXE core *.core OUT*
