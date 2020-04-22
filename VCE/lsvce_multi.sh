args=( "$@" )

echo "${args[@]}" | cut -d " " -f 4- > nset.in

export OMP_NUM_THREADS=18
echo $OMP_NUM_THREADS
gfortran -O3 test_argument.f90 ./final/operationopenmp.f90-o lsvceiter -fopenmp #-heap-arrays -g -traceback -check all -fp-stack-check
time ./lsvceiter $1 $2 $3

