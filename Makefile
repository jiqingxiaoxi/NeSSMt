#!/usr/bin/env make
all:
	gcc -Wall -c NeSSMt.c
	gcc -L /usr/local/lib NeSSMt.o -lgsl -lgslcblas -lm -o NeSSMt
	rm NeSSMt.o

clean:
	rm -f NeSSMt

