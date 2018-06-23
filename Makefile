FC=ifort
CC=gcc

all: bindings_seq
	python -m run.py

bindings_seq:
	F90=$(FC) CC=$(CC); f2py -c -m mylib src/mylib.f90