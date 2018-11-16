#!/bin/bash
exec > _scheduler-stdout.txt
exec 2> _scheduler-stderr.txt
ulimit -v 57042534


export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

'mpirun' '-np' '12' '/home/mborelli/QE_builds/QE_6.3_backports-20181011_no-openmp_oldxml/bin/pw.x' '-nk' '2' '-in' 'aiida.in'  > 'aiida.out' 
