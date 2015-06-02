#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define K_b 9.36E-18
#define sigma 1
#define epsilon 1
#define E 4
#define tot_molecules 20

typedef struct location {
	double xpos;
	double ypos;
	double zpos;
} location;

typedef struct molecule {
	double xpos;
	double ypos;
	double zpos;
	double pos[3];
} molecule;

typedef struct idealGas {
	double size;
	int n_molecules;
	molecule molecules[tot_molecules];
} idealGas;

idealGas *create_gas(double size, int n_molecules) {
	int i;
	double x, y, z;
	srand(time(NULL));
	idealGas *gas = malloc(sizeof(idealGas));

	for (i = 0; i < n_molecules; i++) {
		x = size * ((double)rand() / RAND_MAX);
		y = size * ((double)rand() / RAND_MAX);
		z = size * ((double)rand() / RAND_MAX);
		molecule mol = {.xpos = x, .ypos = y, .zpos = z, .pos = {x, y, z}};
		gas->molecules[i] = mol;
	}

	gas->n_molecules = n_molecules;

	return (gas);
}

double distance(idealGas *gas, int i, int j) {
	int k;

	double distance = 0;
	double *r_i;
	double *r_j;

	r_i = gas->molecules[i].pos;
	r_j = gas->molecules[j].pos;

	for (k = 0; k < 3; k++) {
		distance += (*(r_i + k) - *(r_j + k)) * (*(r_i + k) - *(r_j + k));
	}
	distance = sqrt(distance);
	return distance;
}

double potEnergy(idealGas *gas, int k) {
	double potEnergy = 0;
	int i;
	double r_ki;
	double A;
	for (i = 0; i < gas->n_molecules; i++) {
		if (i != k) {
			r_ki = distance(gas, i, k);	
			A = sigma / r_ki;
			potEnergy += E * (pow(A, 6) - pow(A, 3));
		}
	}
	return potEnergy;
}

double totEnergy(idealGas *gas) {
	double totEnergy = 0;
	int i;
	for (i = 0; i < gas->n_molecules; i++) {
		totEnergy += potEnergy(gas, i);
	}
	return (totEnergy / 2);
}

double randnum(int min, int max)  {
	return(min + (max - min)*((double)rand() / RAND_MAX));
}

void posChange(idealGas *gas, double temperature, int i, double dr, 
	int *changes) {

	double U_i;
	double U_f;
	double dx, dy, dz;
	double x_i, y_i, z_i;

	dx = randnum(-1.0*dr, dr);
	dy = randnum(-1.0*dr, dr);
	dz = randnum(-1.0*dr, dr);

	U_i = totEnergy(gas);

	x_i = gas->molecules[i].xpos;
	y_i = gas->molecules[i].ypos;
	z_i = gas->molecules[i].zpos;

	gas->molecules[i].xpos = x_i + dx;
	gas->molecules[i].ypos = y_i + dy;
	gas->molecules[i].zpos = z_i + dz;

	U_f = totEnergy(gas);

	if (U_f - U_i < 0) {
		*(changes + i) = *(changes + i) + 1;
	} else {
		double prob = exp((-1.0 * (U_f - U_i) / K_b * temperature));
		if (randnum(0, 1) < prob) {
			*(changes + i) = *(changes + i) + 1;
		} else {
			gas->molecules[i].xpos = x_i;
			gas->molecules[i].ypos = y_i;
			gas->molecules[i].zpos = z_i;
		}
	}
	molecule mol = {.xpos = gas->molecules[i].xpos,
		.ypos = gas->molecules[i].ypos,
		.zpos = gas->molecules[i].zpos,
		.pos = {gas->molecules[i].xpos,
				gas->molecules[i].ypos,
				gas->molecules[i].zpos}};
	/*	
	gas->molecules[i].pos = {gas->molecules[i].xpos,
			gas->molecules[i].ypos,
			gas->molecules[i].zpos};
*/
/*
		gas->molecules[i].pos = {.xpos = gas->molecules[i].xpos,
			.ypos = gas->molecules[i].ypos,
			.zpos = gas->molecules[i].zpos};
			*/
}

int posUpdate(idealGas *gas, double temperature, double dr) {
	int changes[gas->n_molecules];

	int i;
	for (i = 0; i < gas->n_molecules; i++) {
		changes[i] = 0;
	}

	int sum;

	for (i = 0; i < gas->n_molecules; i++) {
		posChange(gas, temperature, i, dr, changes);
	}

	for (i = 0; i < gas->n_molecules; i++) {
		sum += *(changes + i);
	}
	return sum;
}

double simulate(idealGas *gas, double temperature, double dr, int trials) {
	int count = 0;
	int i;

	for (i = 0; i < trials; i++) {
		count += posUpdate(gas, temperature, dr);
	}

	return (double)count / (trials * gas->n_molecules);
}

int main(int argc, char *argv[]) {
	srand(time(NULL));
	idealGas *helium = create_gas(20, tot_molecules); 
	double success = simulate(helium, 300, 1, 1000);
	printf("Success rate: %f\n", success);
	return 0;
}
	
