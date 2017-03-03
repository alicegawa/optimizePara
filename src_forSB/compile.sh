make clean
make
mpicc -O3 -o make_neuro_spawn make_neuro_spawn.c
mpicc -O3 -o test_est_target test_est_target.c
