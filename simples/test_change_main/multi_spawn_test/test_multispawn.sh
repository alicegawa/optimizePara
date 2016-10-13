#!/bin/bash -x
#
#PJM -L "rscgrp=fx-xlarge"
#PJM -L "node=6"
#PJM --mpi "shape=1"
#PJM --mpi "proc=32"
#PJM -L "elapse=00:01:00"
#PJM --stg-transfiles all
#--PJM --mpi "use-rankdir"
#PJM --fs large2

#PJM -s

EXEC_NAME="./multi_spawn_test2"

export OMP_NUM_THREADS=32


LPG="lpgparm -t 4MB -s 4MB -d 4MB -p 4MB"
MPIEXEC="mpiexec -np 1 -mca mpi_print_stats 1"
PROF=""
echo "${PROF} ${MPIEXEC} ${LPG} ${EXEC_NAME}"
time ${PROF} ${MPIEXEC} ${LPG} ${EXEC_NAME}

sync