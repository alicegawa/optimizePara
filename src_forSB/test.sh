make clean
rm estimate_main
./compile.sh
num_of_pop=128
mu=64
num_of_child_procs=2
maxeval=100000
num_of_grandchild_procs=4
exec_prog="./test_est_target"
dim_con_mat=3
con_mat_name="../data/conMat_test.txt";
parameter_filename="../data/params_test.txt";
mpiexec -np 1 ./estimate_main ${num_of_pop} ${mu} ${num_of_child_procs} ${maxeval} ${num_of_grandchild_procs} ${exec_prog} ${dim_con_mat} ${con_mat_name} ${parameter_filename}
python judge_results.py
