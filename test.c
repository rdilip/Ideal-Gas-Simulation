#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
double randdob(int min, int max); 
int main(int argc, char *argv[]) {
	printf("modulus is %f\n", fmod(5, 7));
	printf("modulus is %f\n", fmod(10, 7));
}
double randdob(int min, int max)  {
	return(min + (max - min)*((double)rand() / RAND_MAX));
}
