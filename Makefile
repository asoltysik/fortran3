FC=ifort
CC=gcc

bindings_seq:
	F90=$(FC) CC=$(CC); f2py -c -m mylib src/mylib.F90