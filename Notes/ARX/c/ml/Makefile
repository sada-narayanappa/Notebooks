# Modified for tcp migration project

SHELL = /bin/sh

CC = gcc

#SOURCES = datalink.c flow.c main.c tcpip.c tcp.c util.c sysdep.h tcpmig.h
#CONFIG_HEADER = conf.h
PROGRAMS = testml

DEFS = -DHAVE_CONFIG_H -I.
CPPFLAGS = 
LDFLAGS = 
LIBS = -lpcap 
OBJECTS =  main.o linalg.o
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

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

$(PROGRAMS): $(OBJECTS) $(DEPENDENCIES)
	@rm -f $(PROGRAMS)
	$(LINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)

main.o: main.c
linalg.o: linalg.c linalp.h

clean:
	rm -f *.o $(PROGRAMS) core *.core 
