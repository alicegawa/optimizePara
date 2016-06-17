#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TEXT_BUFFER_SIZE 2048
#define MAX_NUM_PARAM 150
#define MAX_NUM_TARGET 5

int main(int argc, char **argv){
    int dimension = 0;
    char buf[TEXT_BUFFER_SIZE];
    double lowerBounds[MAX_NUM_PARAM];
    double upperBounds[MAX_NUM_PARAM];
    unsigned int flg_log[MAX_NUM_PARAM];
    unsigned char flg_log2[MAX_NUM_PARAM];

    FILE *fp;
    if((fp=fopen("./data/params.txt", "r"))==NULL){
	printf("file open eror\n");
	exit(EXIT_FAILURE);
    }
    while( fgets(buf, TEXT_BUFFER_SIZE, fp) != NULL){
	if(strncmp(buf, "#", 1) == 0){ continue;}
	sscanf(buf, "%*s\t%lf\t%lf\t%*lf\t%d\n", &lowerBounds[dimension], &upperBounds[dimension], &flg_log[dimension]);
	//flg_log[dimension] = (unsigned char)atoi((const char*)&flg_log[dimension]);
	printf("flg_log[%d] = %d\n", dimension, flg_log[dimension]);

	sscanf(buf, "%*s\t%lf\t%lf\t%*lf\t%c\n", &lowerBounds[dimension], &upperBounds[dimension], &flg_log2[dimension]);
	flg_log2[dimension] = (unsigned char)atoi((const char*)&flg_log2[dimension]);
	printf("flg_log2[%d] = %d (%c)\n", dimension, flg_log2[dimension], (char)flg_log2[dimension]);
	
	dimension++;//make the dimension information here
    }
    fclose(fp);
    return 0;
}    
