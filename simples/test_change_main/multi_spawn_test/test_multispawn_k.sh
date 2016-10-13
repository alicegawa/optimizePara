#!/bin/bash -x
#
#PJM -L "rscgrp=small"
#PJM -L "node=6"
#PJM --mpi "shape=1"
#PJM --mpi "proc=8"
#PJM -L "elapse=00:01:00"
#PJM --stg-transfiles all
#PJM --stgin "./multi_spawn_test2 ./"
#PJM --stgin "./multi_spawn_test2 0:../"
#PJM --stgin "./call_add ./"
#PJM --stgin "./call_add 0:../"
#PJM --stgin "./add ./"
#PJM --stgin "./add 0:../"
#--PJM --mpi "use-rankdir"

#PJM -s

. /work/system/Env_base

EXEC_NAME="./multi_spawn_test2"

export OMP_NUM_THREADS=8


LPG="lpgparm -t 4MB -s 4MB -d 4MB -p 4MB"
MPIEXEC="mpiexec -n 1 -mca mpi_print_stats 1"
PROF=""
echo "${PROF} ${MPIEXEC} ${LPG} ${EXEC_NAME}"
time ${PROF} ${MPIEXEC} ${LPG} ${EXEC_NAME}

sync
