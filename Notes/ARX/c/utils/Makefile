# Modified for tcp migration project

SHELL = /bin/sh

CC = g++

#SOURCES = datalink.c flow.c main.c tcpip.c tcp.c util.c sysdep.h tcpmig.h
#CONFIG_HEADER = conf.h
PROGRAMS = INVX.EXE
PROGRAMS1 = SCORE.EXE

DEFS = -DHAVE_CONFIG_H -I.
INCLUDES =  -I".." `pkg-config --cflags eigen3`
CPPFLAGS = 
LDFLAGS = 
LIBS = -lpcap 
OBJECTS  =  linreg.o CSV.o model.o LR.o CSVUtils.o main.o invx.o
OBJECTS1 =  CSV.o score.o model.o ncsv.o
DEPENDENCIES = main.o score.o 
CFLAGS = -g -O2 -Wall -Wno-unused
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

all: $(PROGRAMS) $(PROGRAMS1)
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

main.o: main.c
linreg.o: linreg.c linreg.h
CSV.o: CSV.c CSV.h
ncsv.o: ncsv.c ncsv.h any.h
model.o: model.c linreg.h
LR.o: LR.c LR.h
invx.o: invx.c invx.h
score.o: score.c invx.h
CSVUtils.o: CSVUtils.c csv.h

clean:
	rm -f *.o $(PROGRAMS) core *.core OUT*
