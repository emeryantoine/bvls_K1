export OMP_NUM_THREADS=18
echo $OMP_NUM_THREADS
ifort -O3 lsvce.f90 -o lsvceiter -heap-arrays -g -traceback -check all -fp-stack-check -openmp
time ./lsvceiter $1 $2
