#include "SPH_defs.h"

//#include "../cases/RSPH_dambreak_orig/SPH_simulation.h"
//#include "../cases/RSPH_test_orig/SPH_simulation.h"
//#include "../cases/RSPH_test_mdbc/SPH_simulation.h"
//#include "../cases/RSPH_test/SPH_simulation.h"

#include "../cases/ALESPH_dambreak/SPH_simulation.h"
//#include "../cases/ALESPH_test/SPH_simulation.h"
//#include "../cases/ALESPH_test_orig/SPH_simulation.h"
//#include "../cases/ALESPH_test_mdbc/SPH_simulation.h"

//#include "../cases/simpleFlowSr_bt/SPH_simulation.h"

int main(){

	printf("working... \n");

	SPH_simulation SIMULATION;

	SIMULATION.PREP_SIMULATION_DATA(); //move to simulation
	SIMULATION.INIT();
	SIMULATION.RUN();

	return EXIT_SUCCESS;
}
