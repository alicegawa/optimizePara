#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int i;
    int numprocs, myid, new_myid, new_size;
    char command[] = "./add";
    MPI_Comm comm1, comm2, parentcomm;
    char **option;
    char *mpiargv[] = { "aaa", "bbb", "-ccc", NULL };
    char *test[] = {"", "", "",""};
    char filename[100];
    int data[4];
    int color, key;
    FILE *fp;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Comm_get_parent(&parentcomm);

    option = (char **)malloc(sizeof(char *) * 5);

    for(i=0;i<5;i++){
	option[i] = (char *)malloc(sizeof(char) * 3);
    }
    
    key = 0;
    color = myid%4;

    MPI_Comm_split(MPI_COMM_WORLD, color, key, &comm1);
    MPI_Comm_rank(comm1, &new_myid);
    MPI_Comm_size(comm1, &new_size);

    if(new_myid==0){
	sprintf(filename,"data%d.dat",color+1);
	if((fp=fopen(filename,"r"))==NULL){
	    printf("file open error\n");
	    return -1;
	}
	
	for(i=0;i<4;i++){
	    fscanf(fp,"%d", &data[i]);
	    sprintf(option[i],"%d",data[i]);
	}
	fclose(fp);
	option[4] = NULL;
	/* printf("option is"); */
	/* for(i=0;i<4;i++){ */
	/*     printf("%s\n",option[i]); */
	/* } */
    }

    

    //MPI_Comm_spawn(command, MPI_ARGV_NULL, 4, MPI_INFO_NULL, 0, MPI_COMM_SELF, &comm1, MPI_ERRCODES_IGNORE);
    //MPI_Comm_spawn(command, MPI_ARGV_NULL, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    MPI_Comm_spawn(command, option, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    fflush(stdout);
    MPI_Comm_free(&comm2);

    MPI_Finalize();

    return 0;
}
