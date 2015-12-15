#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int n, myid, numprocs, i;
    int new_myid, new_size;
    MPI_Comm comm1, comm2;
    int color, key;
    FILE *fp;
    char filename[100];
    int data[4];
    int num_memofgroup=1, call_final;
    char command[]="./print";
    char **option;
    int maxprocs=2;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    option = (char **)malloc(sizeof(char*) * 4);
    for(i=0;i<4;i++){
	option[i] = (char *)malloc(sizeof(char) * 10);
    }

    if(argc>1){
	num_memofgroup = atoi(argv[1]);
    }
        
    for(i=0;i<4;i++){
	data[i] = 0;
    }
        
    key = 0;
    color = myid/num_memofgroup;
    call_final = numprocs/num_memofgroup;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm1);
    MPI_Comm_rank(comm1, &new_myid);
    MPI_Comm_size(comm1, &new_size);
    
    if(new_myid==0){
	sprintf(filename,"data%d.dat",color+1);
	
	if((fp=fopen(filename,"r"))==NULL){
	    printf("file cannot be open");
	    return -1;
	}

	for(i=0;i<4;i++){
	    fscanf(fp,"%d",&data[i]);
	    sprintf(option[i],"%d",data[i]);
	    //printf("option[%d] is %s\n",i, option[i]);
	    //printf("root of each group color is %d, data[%d] = %d\n", color, i, data[i]);
	}

	fclose(fp);
    }

    printf("test for comm spawn\n");
    //test for comm_spawn
    printf("my original ID = %d, start the mpi comm spawn section\n",myid);
    MPI_Comm_spawn(command, option, 2, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    fflush(stdout);
    MPI_Comm_free(&comm2);
    
    printf("my original ID = %d, finish the mpi comm spawn section\n",myid);
    /*test for broadcast*/
    /* if(new_myid==1){ */
    /* 	for(i=0;i<4;i++){ */
    /* 	    printf("num:1 of each group color is %d, data[%d] = %d\n", color, i, data[i]); */
    /* 	} */
    /* } */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* printf("\n"); */

    /* MPI_Bcast(data, 4, MPI_INT, new_myid, comm1); */
    
    /* if(new_myid==0){ */
    /* 	for(i=0;i<4;i++){ */
    /* 	    printf("root of each group color is %d, data[%d] = %d\n", color, i, data[i]); */
    /* 	} */
    /* } */
    /* if(new_myid==1){ */
    /* 	for(i=0;i<4;i++){ */
    /* 	    printf("num:1 of each group color is %d, data[%d] = %d\n", color, i, data[i]); */
    /* 	} */
    /* } */

    for(i=0;i<4;i++){
	free(option[i]);
    }
    free(option);
    
    MPI_Finalize();
}
