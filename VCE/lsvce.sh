ifort -O3 lsvce.f90 -o lsvceiter -heap-arrays -g -traceback -check all -fp-stack-check
time ./lsvceiter $1
