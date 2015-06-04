#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <signal.h>

#define K_b 9.36E-18
#define SIGMA 1
#define EPSILON 1
#define E 4
#define TOT_MOLECULES 50
#define SIZE 15

void arrPrint(int *ptr, int len);
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
	molecule molecules[TOT_MOLECULES];
} idealGas;

double randnum(double min, double max)  {
	return(min + (max - min)*((double)rand() / RAND_MAX));
}

void die(idealGas *gas, const char *message) {
	if (errno) {
		perror(message);
	} else {
		printf("ERROR: %s\n", message);
	}
	if (gas) {
		free(gas);
	}
	exit(1);
}



double resize(double p, double size) {
	double temp = p;
	if (temp > size) {
		temp = temp- size;
	} else if (p < 0) {
		temp = size + temp;
	}
	
	return temp;
}

double horDist(double x1, double x2, double l) {
	double x_t = x2;
	if (fabs(x1 - x2) > l / 2) {
		x_t = x_t - l;
	}
	return (x1 - x_t);
}


double distance(idealGas *gas, int i, int j) {
	int k;

	double distance = 0;
	double *r_i;
	double *r_j;

	r_i = gas->molecules[i].pos;
	r_j = gas->molecules[j].pos;

	for (k = 0; k < 3; k++) {
		distance += horDist(*(r_i+k),*(r_j+k), SIZE) * 
			horDist(*(r_i+k),*(r_j+k), SIZE);
	}

	distance = sqrt(distance);
	return distance;
}

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

double potEnergy(idealGas *gas, int k, double size) {
	double potEnergy = 0;
	int i;
	double r_ki;
	double A;
	for (i = 0; i < gas->n_molecules; i++) {
		if (i != k) {
			r_ki = distance(gas, i, k);	
			A = SIGMA / r_ki;
			potEnergy += E * (pow(A, 6) - pow(A, 3));
		}
	}
	return potEnergy;
}

double totEnergy(idealGas *gas) {
	double totEnergy = 0;
	int i;
	for (i = 0; i < gas->n_molecules; i++) {
		totEnergy += potEnergy(gas, i, gas->size);
	}
	return (totEnergy / 2);
}

void saveGas(idealGas *gas, const char *filename) {
	int i;
	FILE *f = fopen(filename, "w");

	for (i = 0; i < gas->n_molecules; i++) {
		fprintf(f, "%f\t%f\t%f\n", gas->molecules[i].xpos,
			gas->molecules[i].ypos, gas->molecules[i].zpos);
	}
}

idealGas *loadGas(const char *filename, int n_molecules) {
	FILE *myfile;
	myfile = fopen(filename, "r");

	double x, y, z;
	int i;
	idealGas *gas = malloc(sizeof(idealGas));

	for (i = 0; i < 3 * n_molecules; i++) {
			fscanf(myfile, "%lf", &x);
			fscanf(myfile, "%lf", &y);
			fscanf(myfile, "%lf", &z);
			molecule mol = {.xpos = x, .ypos = y, .zpos = z, .pos = {x, y, z}};
			gas->molecules[i] = mol;
	}

	gas->n_molecules = n_molecules;
	return gas;
}

void printMolecule(molecule mol) {
	printf("(%.15f, %.15f, %.15f)\n", mol.xpos, mol.ypos, mol.zpos);	
}

void printGas(idealGas *gas) {
	int i;
	for (i = 0; i < gas->n_molecules; i++) {
		printMolecule(gas->molecules[i]);
	}
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

	molecule mol = {.xpos = resize(x_i + dx, SIZE),\
		.ypos = resize(y_i + dy , SIZE),\
		.zpos = resize(z_i + dz , SIZE),\
		.pos = {resize(x_i + dx, 10),
				resize(y_i + dy, SIZE),
				resize(z_i + dz , SIZE)}};

	gas->molecules[i] = mol;

	U_f = totEnergy(gas);

	if (U_f - U_i < 0) {
		*(changes + i) = *(changes + i) + 1;
	} else {
		double prob = exp((-1.0 * (U_f - U_i) / K_b * temperature));
		if (randnum(0, 1) < prob) {
			*(changes + i) = *(changes + i) + 1;
		} else {
			molecule mol1 = {.xpos = x_i, .ypos = y_i, .zpos = z_i, 
				.pos = {x_i, y_i, z_i}};
			gas->molecules[i] = mol1;
		}
	}
}

int posUpdate(idealGas *gas, double temperature, double dr) {
	int changes[gas->n_molecules];

	int i;
	for (i = 0; i < gas->n_molecules; i++) {
		changes[i] = 0;
	}

	int sum = 0;

	for (i = 0; i < gas->n_molecules; i++) {
		posChange(gas, temperature, i, dr, changes);
	}

	for (i = 0; i < gas->n_molecules; i++) {
		sum += changes[i];
	}
	return sum;
}

void simulate(idealGas *gas, double T, int num, int data_int,
	double dr, char mode) {
	double energy;
	int i;
	for (i = 0; i < num; i++) {
		posUpdate(gas, T, dr);		
		if (i % data_int == 0) {
			energy = totEnergy(gas);	
			printf("\n\nWrote %d energy: %f\n", i, energy);
			if (mode == 'd') printGas(gas);
		}
	}
}

int main(int argc, char *argv[]) {
	srand(time(NULL));
	char action = argv[1][0];
	idealGas *helium = NULL;
	const char *config;
	double dr, T;
	int num, data_int;

	switch(action) {
		case 'c':
			// Creates new random idealGas, writes positions to file
			// config.txt
			helium = create_gas(SIZE, TOT_MOLECULES); 
			FILE *f1 = fopen("energy_vals1.txt", "w");

			if (f1 == NULL) {
				die(helium, "Error opening file!\n");
			}

			saveGas(helium, "config.txt");
			break;

		case 's':
			// simulates gas molecules given dr and temperature.
			// Also requires number of runs, interval of data
			// collection, and config file.
			if (argc != 7) {
				die(helium, "Case s:Need dr, temp, num runs,"
					" interval of data collection, and"
					" configuration file\n");
			}

			dr = strtod(argv[2], NULL);
			T = strtod(argv[3], NULL);
			num = atoi(argv[4]);
			data_int = atoi(argv[5]);
			config = argv[6];

			if (num < data_int) {
				die(helium, "ERROR: Data collection interval cannot"
					"be greater than total runs\n");
			}

			helium = loadGas(config, TOT_MOLECULES);
			saveGas(helium, "config.txt");
			simulate(helium, T, num, data_int, dr, 'c');	
			saveGas(helium, "final_config.txt");	
			break;

		case 'd':

			// Essentially the same as case s, but also prints
			// out all the gas data each dat_int runs as well
			if (argc != 7) {
				die(helium, "Case s:Need dr, temp, num runs,"
					" interval of data collection, and"
					" configuration file\n");
			}

			dr = strtod(argv[2], NULL);
			T = strtod(argv[3], NULL);
			num = atoi(argv[4]);
			data_int = atoi(argv[5]);
			config = argv[6];

			if (num < data_int) {
				die(helium, "ERROR: Data collection interval cannot"
					"be greater than total runs\n");
			}

			helium = loadGas(config, TOT_MOLECULES);
			saveGas(helium, "config.txt");
			simulate(helium, T, num, data_int, dr, 'd');	
			saveGas(helium, "final_config.txt");	
			break;

		default:
			die(helium, "Defualt: Need dr, temp, num runs, interval of "
				"data collection, and configuration file\n");
		}
			
	free(helium);

	return 0;
}
