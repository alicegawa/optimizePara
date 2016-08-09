#!/bin/bash -x
#
#PJM --rsc-list "rscgrp=large"
#PJM --rsc-list "node=9216"
#PJM --rsc-list "elapse=02:00:00"
#PJM --mpi "shape=1024"
#PJM --mpi "proc=8192"
#PJM --stg-transfiles all
#PJM --stgin "rank=* ./estimate_main %r:./"
#PJM --stgin "rank=0 ./estimate_main %r:../"
#PJM --stgin "rank=* ./*.err %r:./"
#PJM --stgin "rank=* ./*.par %r:./"
#PJM --stgin "rank=* ../hocfile/* %r:./"
#PJM --stgin "rank=0 ../hocfile/* %r:../"
#PJM --stgin "rank=* ../data/* %r:./data/"
#PJM --stgin "rank=0 ../data/* %r:../data/"
#PJM --stgin "rank=* ../data/rev_potential.dat %r:./"
#PJM --stgin "rank=0 ../data/rev_potential.dat %r:../"
#PJM --stgin "rank=* ../data/params.txt %r:./"
#PJM --stgin "rank=0 ../data/params.txt %r:../"


#PJM --stgin "rank=* /home/hp120263/k01792/neuron_kplus/specials/sparc64/special %r:./"
#PJM --stgin "rank=0 /home/hp120263/k01792/neuron_kplus/specials/sparc64/special %r:../"
#PJM --stgin "rank=* /home/hp120263/k01792/neuron_kplus/stgin/* %r:./"
#PJM --stgin "rank=0 /home/hp120263/k01792/neuron_kplus/stgin/* %r:../"

#PJM --mpi "use-rankdir"
#PJM --stgout "rank=* %r:./filib_rtinfo* ./"
#PJM --stgout "rank=0 %r:../flib_rtinfo* ./"

#PJM -s
#PJM -m "e"
#PJM --mail-list "fukuda@g.brain.imi.i.u-tokyo.ac.jp"
#

. /work/system/Env_base

#export PARALLEL=256
#export OMP_NUM_THREADS=$PARALLEL

export FLIB_RTINFO_PA=Statistics
export FLIB_RTINFO_PROCESS=txt

ls
printenv

MPIEXEC="mpiexec -mca mpi_print_stats 1"
NPROC="-n 8192"
PROF=""
EXEC_FILE="./estimate_main"
echo "${PROF} ${MPIEXEC} ${NPROC} ${EXEC_FILE}"
time  ${PROF} ${MPIEXEC} ${NPROC} ${EXEC_FILE}

sync