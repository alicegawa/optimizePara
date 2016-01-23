#!/bin/bash -x
#
#PJM --rsc-list "node=9"
#PJM --rsc-list "elapse=00:01:00"
#PJM --mpi "shape=1"
#PJM --mpi "proc=8"
#PJM --stg-transfiles all
#PJM --stgin "./print_color_id ./"
#PJM --stgin "./spawn_print_k ./"
#PJM -s
#

. /work/system/Env_base

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL


MPIEXEC="mpiexec -mca mpi_print_stats 1"
NPROC="-n 8"
PROF=""
EXEC_FILE="./spawn_print_k"
echo "${PROF} ${MPIEXEC} ${NPROC} ${EXEC_FILE}"
time ${PROF} ${MPIEXEC} ${NPROC} ${EXEC_FILE}

sync