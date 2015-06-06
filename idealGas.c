#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>
#include <signal.h>

#define K_b 0.0026994
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
	double x12 = x2 - x1;
	return (x12 - l * round(x12/l));
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

idealGas *create_gas(double size) {
	int i;
	double x, y, z;
	srand(time(NULL));
	idealGas *gas = malloc(sizeof(idealGas));

	for (i = 0; i < TOT_MOLECULES; i++) {
		x = size * ((double)rand() / RAND_MAX);
		y = size * ((double)rand() / RAND_MAX);
		z = size * ((double)rand() / RAND_MAX);
		molecule mol = {.xpos = x, .ypos = y, .zpos = z, .pos = {x, y, z}};
		gas->molecules[i] = mol;
	}
	gas->n_molecules = TOT_MOLECULES;

	return (gas);
}

double potEnergy(idealGas *gas, int k, double size) {
	double potEnergy = 0;
	int i;
	double r_ki;
	double A;
	for (i = 0; i < TOT_MOLECULES; i++) {
		if (i != k) {
			r_ki = distance(gas, i, k);	
			A = SIGMA / r_ki;
			potEnergy += E * (pow(A, 6) - 2 * pow(A, 3));
		}
	}
	return potEnergy;
}

double totEnergy(idealGas *gas) {
	double totEnergy = 0;
	int i;
	for (i = 0; i < TOT_MOLECULES; i++) {
		totEnergy += potEnergy(gas, i, gas->size);
	}
	return (totEnergy / 2);
}

void saveGas(idealGas *gas, const char *filename) {
	int i;
	FILE *f = fopen(filename, "w");

	for (i = 0; i < TOT_MOLECULES; i++) {
		fprintf(f, "%f\t%f\t%f\n", gas->molecules[i].xpos,
			gas->molecules[i].ypos, gas->molecules[i].zpos);
	}
	fflush(f);
	fclose(f);
}

idealGas *loadGas(const char *filename) {
	FILE *myfile;
	myfile = fopen(filename, "r");

	double x, y, z;
	int i;
	idealGas *gas = malloc(sizeof(idealGas));

	for (i = 0; i < 3 * TOT_MOLECULES; i++) {
			fscanf(myfile, "%lf", &x);
			fscanf(myfile, "%lf", &y);
			fscanf(myfile, "%lf", &z);
			molecule mol = {.xpos = x, .ypos = y, .zpos = z, .pos = {x, y, z}};
			gas->molecules[i] = mol;
	}

	gas->n_molecules = TOT_MOLECULES;
	fclose(myfile);
	return gas;
}

void printMolecule(molecule mol) {
	printf("(%.4f, %.4f, %.4f)\n", mol.xpos, mol.ypos, mol.zpos);	
}

void printGas(idealGas *gas) {
	int i;
	for (i = 0; i < TOT_MOLECULES; i++) {
		printMolecule(gas->molecules[i]);
	}
}

int posChange(idealGas *gas, double temperature, int i, double dr) {
	double U_i;
	double U_f;
	double dx, dy, dz;
	double x_i, y_i, z_i;
	int res = 0;

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
		res = 1;
	} else {
		double prob = exp((-1.0 * (U_f - U_i) / K_b * temperature));
		if (randnum(0, 1) < prob) {
			res = 1;
		} else {
			molecule mol1 = {.xpos = x_i, .ypos = y_i, .zpos = z_i, 
				.pos = {x_i, y_i, z_i}};
			gas->molecules[i] = mol1;
		}
	}
	return res;
}

int posUpdate(idealGas *gas, double temperature, double dr) {

	int i;
	int sum = 0;

	for (i = 0; i < TOT_MOLECULES; i++) {
		sum += posChange(gas, temperature, i, dr);
	}
	return sum;
}

void simulate(idealGas *gas, double T, int num, int data_int,
	double dr, char mode) {
	double energy;
	FILE *f = fopen("energy.txt", "w");
	int i;
	for (i = 0; i < num; i++) {
		posUpdate(gas, T, dr);		
		if (i % data_int == 0) {
			if (mode == 'd') printGas(gas);
			energy = totEnergy(gas);	
			printf("Wrote %d energy: %f\n", i, energy);
			fprintf(f, "%f\n", energy);
			fflush(f);
			saveGas(gas, "final_config.txt");
		}
	}
	fclose(f);
}

double findAcceptanceRate(const char *config, double T, double dr, int num) {
	int count = 0;
	int i;
	idealGas *gas = loadGas(config);
	for (i = 0; i < num; i++) {
		count += posUpdate(gas, T, dr);
	}

	printf("%d\n", count);
	return ((double)count / (double)(num * TOT_MOLECULES));
}

void equilibriate(const char *config, double Ti, double Tf, double dr) {
	double T2 = 0.66 * (Ti - Tf) + Tf;
	double T3 = 0.33 * (Ti - Tf) + Tf;
	int trials = 3000;
	idealGas *gas = loadGas(config);	
	int i;
	for (i = 0; i < trials; i++) {
		posUpdate(gas, Ti, dr);
	}
	printf("Ran with temperature %.5f K\n", Ti);
	saveGas(gas, "final_config.txt");
	for (i = 0; i < trials; i++) {
		posUpdate(gas, T2, dr);
	}
	printf("Ran with temperature %.5f K\n", T2);
	saveGas(gas, "final_config.txt");
	for (i = 0; i < trials; i++) {
		posUpdate(gas, T3, dr);
	}
	printf("Ran with temperature %.5f K\n", T3);

	saveGas(gas, "final_config.txt");
	for (i = 0; i < trials; i++) {
		posUpdate(gas, Tf, dr);
	}
	printf("Ran with temperature %.5f K\n", Tf);

	saveGas(gas, "final_config.txt");

}

int main(int argc, char *argv[]) {
	srand(time(NULL));
	char action = argv[1][0];
	idealGas *helium = NULL;
	const char *config;
	double dr, T, Ti, Tf;
	int num, data_int;

	switch(action) {
		case 'c':
			// Creates new random idealGas, writes positions to file
			// config.txt
			helium = create_gas(SIZE); 
			FILE *f1 = fopen("energy_vals1.txt", "w");

			if (f1 == NULL) {
				die(helium, "Error opening file!\n");
			}

			saveGas(helium, "config.txt");
			fclose(f1);
			break;
		case 'e':
			// Equilibriates gas
			if (argc != 6) {die(helium, "Need dr, Ti, Tf, and configuration"
				" file to equilibriate");}

			dr = strtod(argv[2], NULL);
			Ti = strtod(argv[3], NULL);
			Tf = strtod(argv[4], NULL);
			config = argv[5];

			equilibriate(config, Ti, Tf, dr);
			break;
		case 'a':
			// Finds acceptance ratio. Should be preceded by equilibriation
			if (argc != 6) {die(helium, "Need temperature, dr, number of "
				"trials, and configuration file to find acceptance rate");}

			T = strtod(argv[2], NULL);
			dr = strtod(argv[3], NULL);
			num = atoi(argv[4]);
			config = argv[5];
			printf("Acceptance rate: %f%%\n", 100 * \
				findAcceptanceRate(config, T,dr, num));
			break;
		case 's':
			// simulates gas molecules given dr and temperature.
			// Also requires number of runs, interval of data
			// collection, and config file.
			if (argc != 7) {
				die(helium, "Need dr, temp, num runs,"
					" interval of data collection, and"
					" configuration file\n");
			}

			dr = strtod(argv[2], NULL);
			T = strtod(argv[3], NULL);
			num = atoi(argv[4]);
			data_int = atoi(argv[5]);
			config = argv[6];

			if (num < data_int) {
				die(helium, "Data collection interval cannot"
					"be greater than total runs\n");
			}

			helium = loadGas(config);
			saveGas(helium, "config.txt");
			simulate(helium, T, num, data_int, dr, 'c');	
			saveGas(helium, "final_config.txt");	
			break;

		case 'd':

			// Essentially the same as case s, but also prints
			// out all the gas data each dat_int runs as well
			if (argc != 7) {
				die(helium, "Need dr, temp, num runs,"
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

			helium = loadGas(config);
			saveGas(helium, "config.txt");
			simulate(helium, T, num, data_int, dr, 'd');	
			saveGas(helium, "final_config.txt");	
			break;

		default:
			die(helium, "Need dr, temp, num runs, interval of "
				"data collection, and configuration file\n");
		}
			
	free(helium);

	return 0;
}
