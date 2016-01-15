#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv){
    int num_cell;
    float **weight, **delay, **rev_e;
    int i, j;
    FILE *fp1, *fp2, *fp3;

    srand((unsigned)time(NULL));
    
    if(argc>1){
	num_cell = atoi(argv[1]);
    }else{
	num_cell = 86;
    }

    weight = (float **)malloc(sizeof(float *) * num_cell);
    delay = (float **)malloc(sizeof(float *) * num_cell);
    rev_e = (float **)malloc(sizeof(float *) * num_cell);
    for(i=0;i<num_cell;i++){ 
	weight[i] = (float *)malloc(sizeof(float) * num_cell);
	delay[i] = (float *)malloc(sizeof(float) * num_cell);
	rev_e[i] = (float *)malloc(sizeof(float) * num_cell);
    }

    for(i=0;i<num_cell;i++){
	for(j=0;j<num_cell;j++){
	    weight[i][j] = 0;
	    rev_e[i][j] = 10;
	    delay[i][j] = (rand()%1000) * 0.02 + 1;
	}
    }

    if((fp1=fopen("weight.dat","w"))==NULL || (fp2=fopen("delay.dat","w"))==NULL || (fp3=fopen("rev_potential.dat","w"))==NULL){
	printf("file open error\n");
	return -1;
    }

    for(i=0;i<num_cell;i++){
	for(j=0;j<num_cell;j++){
	    fprintf(fp1,"%7.4f\t",weight[i][j]);
	    fprintf(fp2,"%7.4f\t",delay[i][j]);
	    fprintf(fp3,"%7.4f\t",rev_e[i][j]);
	}
	fprintf(fp1,"\n");
	fprintf(fp2,"\n");
	fprintf(fp3,"\n");
    }

    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    return 0;
}
