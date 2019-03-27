# Modified for tcp migration project

SHELL = /bin/sh

OS=$(shell uname -s)
CC = g++

#SOURCES = datalink.c flow.c main.c tcpip.c tcp.c util.c sysdep.h tcpmig.h
#CONFIG_HEADER = conf.h
PROGRAMS1 = SCORE.EXE_$(OS)

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
#AM_CFLAGS=-static-libstdc++ -static-libgcc
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


main.o: main.c
CSV.o: CSV.c CSV.h
ncsv.o: ncsv.c ncsv.h any.h
score.o: score.c any.h
CSVUtils.o: CSVUtils.c CSV.h


clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS1) INV.EXE SCORE.EXE core *.core OUT*
