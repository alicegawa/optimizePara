#include<stdio.h>

int main(int argc, char **argv){
    int num;
    if(argc<2){
	printf("suffix is not enough\n");
	return -1;
    }

    num = atoi(argv[1]);
    printf("input number is %d\n", num);
    return 0;
}
