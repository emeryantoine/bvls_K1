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
src = lsvce.f90 ../operationopenmp.f90
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

help:
	@echo "file : $(file)" 
	@echo "size : $(size)"
	@echo "nbr_threads : $(nbr_threads)"
	@echo "CC : $(CC)"
	@echo "flags : $(flags)"
	@echo "exec : $(exec)"
	@echo "src : $(src)"
	@echo "make, make run, make time, make help, make clean"

clean:
	rm -rf $(exec)
