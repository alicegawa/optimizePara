#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

void read_answers(double *ref, int main_myid, int gene_dim, int true_dim, int eval_dim){
  int i;
  FILE *fp;
  double *ref_full=NULL;
  
  if((fp=fopen("answer.txt", "r"))==NULL){
    printf("file open error occurs in read_answers\n");
  }
  ref_full = (double *)malloc(sizeof(double) * gene_dim);
  for(i=0;i<gene_dim;++i){
    if(i<true_dim){
      fscanf(fp, "%lf", &ref_full[i]);
    }else{
      ref_full[i] = 0;
    }
  }
  for(i=0;i<eval_dim;++i){
    ref[i] = ref_full[eval_dim * main_myid + i];
  }
  fclose(fp);
  free(ref_full);
}

double evaluation(double *gene_info, double *ref, int gene_dim){
  int i;
  double score=0.0;
  for(i=0;i<gene_dim;++i){
    score += (gene_info[i] - ref[i]) * (gene_info[i] - ref[i]);
  }
  return score;
}

int main(int argc, char **argv){
  int i;
  MPI_Comm parentcomm, spawn_parent_comm;
  int main_size, main_myid;
  int spawn_parent_size, spawn_parent_rank;

  int num_of_my_pop;
  int gene_dim = 36;
  int true_dim;
  int eval_dim_per_proc;
  double *gene_info=NULL, *ref=NULL;
  double *gene_info_full;/*defined only in parent procs*/
  double *scores=NULL, *scores_reduce=NULL;
  double *scores_full=NULL;
  
  double flg_termination=0.0;
  double fbest=10000;
  int fbest_idx;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &main_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &main_myid);

  /* make new intracommunicator for communicating with the parent process*/
  MPI_Comm_get_parent(&parentcomm);
  MPI_Intercomm_merge(parentcomm, 1, &spawn_parent_comm);
  MPI_Comm_size(spawn_parent_comm, &spawn_parent_size);
  MPI_Comm_rank(spawn_parent_comm, &spawn_parent_rank);

  /*receive info from argv?*/
  if(argc<2){
    printf("please type the problem setting properly\n");
    printf("set \"num_of_my_pop\" and \"gene_dim\" default value\n");
    num_of_my_pop = 1;
    gene_dim = 36;
    true_dim = 36;
  }else{
    num_of_my_pop = atoi(argv[1]);
    gene_dim = atoi(argv[2]);
    true_dim = atoi(argv[3]);
  }
  //printf("info: num_of_my_pop=%d, gene_dim=%d, true_dim = %d, main_size=%d\n", num_of_my_pop, gene_dim, true_dim, main_size);
  eval_dim_per_proc = gene_dim / main_size;
  //printf("eval_dim_per_proc = %d\n", eval_dim_per_proc);

  /*alloc the memory for estimation*/
  gene_info = (double *)malloc(sizeof(double) * eval_dim_per_proc * num_of_my_pop);
  if(gene_info==NULL){
    printf("memory allocation error occurs @{gene_info}\n");
  }
  ref = (double *)malloc(sizeof(double) * eval_dim_per_proc);
  if(ref==NULL){
    printf("memory allocation error occurs @{ref}\n");
  }
  scores = (double *)malloc(sizeof(double) * num_of_my_pop);
  if(scores==NULL){
    printf("memory allocation error occurs @{scores}\n");
  }
  scores_reduce = (double *)malloc(sizeof(double) * num_of_my_pop);
  if(scores_reduce==NULL){
    printf("memory allocation error occurs @{scores_reduce}\n");
  }
  scores_full = (double *)calloc(num_of_my_pop * spawn_parent_size, sizeof(double));
  if(scores_full==NULL){
    printf("memory allocation error occurs @{scores_full}\n");
  }

  /* read setting file*/
  read_answers(ref, main_myid, gene_dim, true_dim, eval_dim_per_proc);

  /* for(i=0;i<gene_dim;++i){ */
  /*   printf("ref[%d] = %lf\t", i, ref[i]); */
  /* } */
  /* printf("\n"); */
  /*get the gene information from parent, evaluate the fitness, and send them to parent*/
  while(1){
    fbest = 10000;
    MPI_Scatter(gene_info_full, eval_dim_per_proc * num_of_my_pop, MPI_DOUBLE, gene_info, eval_dim_per_proc * num_of_my_pop, MPI_DOUBLE, 0, spawn_parent_comm);
    for(i=0; i<num_of_my_pop; ++i){
      scores[i] = evaluation(&gene_info[i*eval_dim_per_proc], ref, eval_dim_per_proc);
    }
    MPI_Allreduce(scores, scores_reduce, num_of_my_pop, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    /* for(i=0;i<num_of_my_pop;++i){ */
    /*   if(fbest > scores_reduce[i]){ */
    /* 	fbest = scores_reduce[i]; */
    /* 	fbest_idx = i; */
    /*   } */
    /* } */
    /* if(fbest < 1){ */
    /*   printf("fbest=%lf\n", fbest); */
    /*   for(i=0;i<eval_dim_per_proc;++i){ */
    /* 	printf("@%d proc, gene_value = %lf(fbest=%lf)\n", main_myid, gene_info[fbest_idx * eval_dim_per_proc+i], fbest); */
    /*   } */
    /* } */

    MPI_Gather(scores_reduce, num_of_my_pop, MPI_DOUBLE, scores_full, num_of_my_pop, MPI_DOUBLE, 0, spawn_parent_comm);
    MPI_Bcast(&flg_termination, 1, MPI_DOUBLE, 0, spawn_parent_comm);
    //printf("flg_termination@{test_est_target} is %lf\n", flg_termination);
    if((int)flg_termination){
      break;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  printf("end of loop\n");

  free(gene_info); gene_info=NULL;
  free(ref); ref=NULL;
  free(scores); scores=NULL;
  free(scores_full); scores_full=NULL;

  MPI_Finalize();
  return 0;
}
