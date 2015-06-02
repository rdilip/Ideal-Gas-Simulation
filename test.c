#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
double randdob(int min, int max); 
int main(int argc, char *argv[]) {
	double x, y, z;
	srand(time(NULL));
	int i;
	for (i = 0; i < 10; i++) {
		printf("%f\n", randdob(1, 4));
	}
}
double randdob(int min, int max)  {
	return(min + (max - min)*((double)rand() / RAND_MAX));
}
