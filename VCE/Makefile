SHELL := /bin/bash

#Variables a modifier si besoin
###########################################
file ?= ../../../transfert/RA.out.mexmro.12
size ?= 12
set ?= 2
setsize ?= 1 2 3 4 5
nbr_threads ?= 18
schedule ?= dynamic,5
###########################################

#Variables du makefile
###########################################
CC = ifort
flags = -O3 -openmp
exec = lsvceiter
src = test_arguments.f90 operationopenmp.f90
###########################################

all:$(exec)

$(exec):$(src)
	$(CC) -o $@ $^ $(flags)

run:$(exec)
	export OMP_NUM_THREADS=$(nbr_threads)
	export OMP_SCHEDULE=$(schedule)
	./$(exec) $(size) $(set) $(file) $(setsize)

time:$(exec)
	export OMP_NUM_THREADS=$(nbr_threads)
	export OMP_SCHEDULE=$(schedule)
	time ./$(exec) $(size) $(set) $(file) $(setsize)

help:
	@echo "file : $(file)" 
	@echo "size : $(size)"
	@echo "nbr_threads : $(nbr_threads)"
	@echo "CC : $(CC)"
	@echo "flags : $(flags)"
	@echo "exec : $(exec)"
	@echo "src : $(src)"
	@echo "set : $(set)"
	@echo "setsize : $(setsize)"
	@echo "make, make run, make time, make help, make clean"

clean:
	rm -rf $(exec)
