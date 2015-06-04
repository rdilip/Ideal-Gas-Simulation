#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

double distance(double *pos1,double *pos2, int i, int j) {
	int k;

	double distance = 0;
	double *r_i;
	double *r_j;

	r_i = pos1;
	r_j = pos2;

	for (k = 0; k < 3; k++) {
		distance += (*(r_i + k) - *(r_j + k)) * (*(r_i + k) - *(r_j + k));
	}

	distance = sqrt(distance);
	return distance;
}
double resize(double p, double size) {
	if (p > size) {
		p = p - size;
	} else if (p < 0) {
		p = size + p;
	}

	return p;
}
double randdob(int min, int max); 
int main(int argc, char *argv[]) {
	foo(7);
}

double randdob(int min, int max)  {
	return(min + (max - min)*((double)rand() / RAND_MAX));
}

