SUBDIRS=ml ml1

all: ins
	cd ml;make;
	cd utils;make;

main: helloworld.o
helloworld.o:
	gcc helloworld.c

PWD=`pwd`

ins:
	echo "****"
	#mkdir -p ~/tcpmig;D=`pwd`;ln -s $(PWD)/scripts ~/tcpmig/scripts
