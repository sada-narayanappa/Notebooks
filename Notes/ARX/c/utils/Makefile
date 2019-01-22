# Modified for tcp migration project

SHELL = /bin/sh

CC = g++

#SOURCES = datalink.c flow.c main.c tcpip.c tcp.c util.c sysdep.h tcpmig.h
#CONFIG_HEADER = conf.h
PROGRAMS = INVX.EXE

DEFS = -DHAVE_CONFIG_H -I.
INCLUDES =  -I".." `pkg-config --cflags eigen3`
CPPFLAGS = 
LDFLAGS = 
LIBS = -lpcap 
OBJECTS =  main.o linreg.o CSV.o model.o LR.o
DEPENDENCIES = 
CFLAGS = -g -O2 -Wall -Wno-unused
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

all: $(PROGRAMS)
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

$(PROGRAMS): $(OBJECTS) $(DEPENDENCIES)
	@rm -f $(PROGRAMS)
	$(LINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)

main.o: main.c
linreg.o: linreg.c linreg.h
CSV.o: CSV.c CSV.h
model.o: model.c linreg.h
LR.o: LR.c LR.h

clean:
	rm -f *.o $(PROGRAMS) core *.core OUT*
