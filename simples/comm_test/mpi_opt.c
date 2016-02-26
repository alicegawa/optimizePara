#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>

int main(int argc, char **argv){
    int n, myid, numprocs, i;
    int new_myid, new_size;
    MPI_Comm comm1, comm2, parentcomm, nrn_comm;
    int color, key;
    FILE *fp;
    char filename[100];
    int data[4];
    int num_memofgroup=1, call_final;
    char command[]="./add";
    char specials[]="special";
    char **option;
    char *neuron_argv[] = {"-mpi", "-nobanner", "para_test.hoc", NULL};
    char **neuron_argv2;
    char option_mpi[] ="-mpi", HOCFILE[] = "para_test.hoc", option_nobanner[]="-nobanner";
    int maxprocs=2;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    //MPI_Comm_get_parent(&parentcomm);

    option = (char **)malloc(sizeof(char*) * 4);
    for(i=0;i<4;i++){
	option[i] = (char *)malloc(sizeof(char) * 10);
    }    

    neuron_argv2 = (char **)malloc(sizeof(char *) * 6);
    for(i=0;i<6;i++){
	neuron_argv2[i] = (char *)malloc(sizeof(char) * 10);
	if(i==0){
	    sprintf(neuron_argv2[i],"%s",option_mpi);
	}else if(i==1){
	    sprintf(neuron_argv2[i],"%s",option_nobanner);
	}else if(i==2){
	    sprintf(neuron_argv2[i],"-c");
	}else if(i==4){
	    sprintf(neuron_argv2[i],"%s",HOCFILE);
	}else{
	}
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

    sprintf(neuron_argv2[3],"COLOR=%d",color);
    neuron_argv2[5] = NULL;
    
    /* if(new_myid==0){ */
    /* 	sprintf(filename,"data%d.dat",color+1); */
	
    /* 	if((fp=fopen(filename,"r"))==NULL){ */
    /* 	    printf("file cannot be open"); */
    /* 	    return -1; */
    /* 	} */

    /* 	for(i=0;i<4;i++){ */
    /* 	    fscanf(fp,"%d",&data[i]); */
    /* 	    sprintf(option[i],"%d",data[i]); */
    /* 	    //printf("option[%d] is %s\n",i, option[i]); */
    /* 	    //printf("root of each group color is %d, data[%d] = %d\n", color, i, data[i]); */
    /* 	} */

    /* 	fclose(fp); */
    /* } */

    printf("test for comm spawn\n");
    //test for comm_spawn
    //printf("my original ID = %d, start the mpi comm spawn section\n",myid);
    //MPI_Comm_spawn(command, option, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    //printf("finish mpi spawn test for \"add\"\n");
    //MPI_Comm_spawn(specials, neuron_argv, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    MPI_Comm_spawn(specials, neuron_argv2, 4, MPI_INFO_NULL, 0, comm1, &comm2, MPI_ERRCODES_IGNORE);
    MPI_Intercomm_merge(comm2,0, &nrn_comm);
    printf("finish comm spawn\n");
    fflush(stdout);
    
    for(i=0;i<4;i++){
	free(option[i]);
    }
    free(option);

    MPI_Finalize();
}
