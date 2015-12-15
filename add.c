#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char** argv){
    int num[4];
    int sum=0;
    int i;
    
    if(argc<5){
	printf("suffix is not enough\n");
	return -1;
    }

    for(i=0;i<4;i++){
	num[i] = atoi(argv[i+1]);
	sum += num[i];
    }

    printf("sum = %d\n",sum);

    return 0;
}
    
