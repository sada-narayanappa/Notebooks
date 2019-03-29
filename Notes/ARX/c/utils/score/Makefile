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
PROGRAMS = SCORE.EXE_$(FLAVOR)
DEFS = -DHAVE_CONFIG_H -I.
INCLUDES =  -I".."
CPPFLAGS = 
LDFLAGS = 
#LIBS = -lpcap  -lpthread
OBJECTS =  CSV.o score.o ncsv.o
DEPENDENCIES = 
CFLAGS = -g -O2 -Wall -Wno-unused
COMPILE = $(CC) -static $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
CCLD = $(CC)
ifeq ($(FLAVOR),alpine)
	AM_CFLAGS=-static-libstdc++ -static-libgcc
endif
LINK = $(CCLD)  $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@

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
	$(LINK) $(OBJECTS) $(LIBS)

lib: $(OBJECTS) $(DEPENDENCIES)
	$(CCLD) -shared -o libscore.so $(OBJECTS) 
    
#$(CCLD) -L. -Wall -o test main.c -lscore

CSV.o: CSV.c CSV.h marray.h
ncsv.o: ncsv.c ncsv.h any.h marray.h
score.o: score.c any.h marray.h

clean:
	rm -f *.o $(PROGRAMS) SCORE.EXE* libscore* core *.core OUT*
