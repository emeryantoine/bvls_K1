SHELL := /bin/bash

#Variables a modifier si besoin
###########################################
file ?= ../../../transfert/RA.out.mexmro.12
size ?= 12
nbr_threads ?= 18
###########################################

#Variables du makefile
###########################################
CC = ifort
flags = -O3 -heap-arrays -g -traceback -check all -fp-stack-check -openmp
exec = lsvceiter
src = lsvce.f90
###########################################

all:$(exec)

$(exec):$(src)
	$(CC) -o $@ $^ $(flags)

run:$(exec)
	export OMP_NUM_THREADS=$(nbr_threads)
	./$(exec) $(size) $(file)

time:$(exec)
	export OMP_NUM_THREADS=$(nbr_threads)
	time ./$(exec) $(size) $(file)

clean:
	rm -rf $(exec)
