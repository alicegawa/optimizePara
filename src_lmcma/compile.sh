mpicc -O3 -c boundary_transformation.c
mpicc -O3 -c my_boundary_transformation.c
mpicc -O3 -c estimate_main.c -lm
mpicc -O3 -o estimate_main boundary_transformation.o my_boundary_transformation.o estimate_main.o -lm
mpicc -O3 -o make_neuro_spawn make_neuro_spawn.c
