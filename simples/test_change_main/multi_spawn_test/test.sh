mpicc -O3 -o call_add call_add.c
mpicc -O3 -o multi_spawn_test2 multi_spawn_test2.c
mpiexec -np 1 ./multi_spawn_test2
