#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

void sig_handler(int signo) {
	if (signo == SIGINT) {
		printf("Received SIGINT\n");
		exit(1);
	}
}

int main(int argc, char *argv[]) {
	if (signal(SIGINT, sig_handler) == SIG_ERR) {
		printf("\ncan't catch yo SIGINT\n");
	}
	while (1) {
		sleep(1);
	}
	return 1;
}
