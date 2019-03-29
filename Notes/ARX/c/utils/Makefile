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


PROGRAMS = INVX.EXE_$(FLAVOR)

CC = g++
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

ifeq ($(FLAVOR),alpine)                    
        AM_CFLAGS=-static-libstdc++ -static-libgcc     
	INCLUDES =  -I".." -I"/usr/include/eigen3" 
endif                   

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

$(PROGRAMS): $(OBJECTS) CSVUtils.o $(DEPENDENCIES)
	@rm -f $(PROGRAMS)
	$(LINK) $(LDFLAGS) $(OBJECTS) $(LDADD) $(LIBS)

$(PROGRAMS1): $(OBJECTS1) $(DEPENDENCIES)
	@rm -f $(PROGRAMS1)
	$(LINK) $(LDFLAGS) $(OBJECTS1) $(LDADD) $(LIBS)

main.o: main.c
CSV.o: CSV.c CSV.h
LR.o: LR.c LR.h
invx.o: invx.c invx.h any.h
CSVUtils.o: CSVUtils.c CSV.h


clean:
	rm -f *.o $(PROGRAMS) $(PROGRAMS1) INVX.EXE* core *.core OUT*
