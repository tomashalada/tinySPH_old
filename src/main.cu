#include "SPH_defs.h"
#include "../cases/flowDrop_bt/SPH_simulation.h"

int main(){

	printf("working... \n");

	SPH_simulation SIMULATION;

	SIMULATION.PREP_SIMULATION_DATA(); //move to simulation
	SIMULATION.INIT();
	SIMULATION.RUN();

	return EXIT_SUCCESS;
}
