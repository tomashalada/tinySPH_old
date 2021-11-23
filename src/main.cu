#include "SPH_defs.h"
#include "/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/ALESPH_dambreakLong/MUSCL.h"

//===================================================================//
//===========================TEST=CASES==============================//
//===================================================================//
// ========= DamBreak mDBC ======== //
//#include "../cases/damBreakLong_mDBC/SPH_simulationSymplectic.h"
//#include "../cases/damBreakLong_mDBC/SPH_simulationVerlet.h"
//#include "../cases/damBreakLong_mDBC/SPH_simulationLeapFrog.h"

// ========= DamBreak DBC ======== //
//#include "../cases/damBreakLong_DBC/SPH_simulationSymplectic.h"
//#include "../cases/damBreakLong_DBC/SPH_simulationVerlet.h"
//#include "../cases/damBreakLong_DBC/SPH_simulationLeapFrog.h"

// ========= DamBreak BT ======== //
//#include "../cases/damBreakLong_BT/SPH_simulationSymplectic.h"
//#include "../cases/damBreakLong_BT/SPH_simulationVerlet.h" /* aktualne pada */
//#include "../cases/damBreakLong_BT/SPH_simulationLeapFrog.h" /* aktualne pada */

//===================================================================//
// ====== StillWatter mDBC ======= //
//#include "../cases/stillWater_mDBC/SPH_simulationSymplectic.h"
//#include "../cases/stillWater_mDBC/SPH_simulationVerlet.h"
//#include "../cases/stillWater_mDBC/SPH_simulationLeapFrog.h"

// ====== StillWatter DBC ======= //
//#include "../cases/stillWater_DBC/SPH_simulationSymplectic.h"
//#include "../cases/stillWater_DBC/SPH_simulationVerlet.h"
//#include "../cases/stillWater_DBC/SPH_simulationLeapFrog.h"

// ====== StillWatter BT ======= //
//#include "../cases/stillWater_BT/SPH_simulationSymplectic.h"
//#include "../cases/stillWater_BT/SPH_simulationVerlet.h"

//===================================================================//
// ======= flowDrop mDBC ======== //
//#include "../cases/flowDrop_mDBC/SPH_simulationSymplectic.h"
//#include "../cases/flowDrop_mDBC/SPH_simulationLeapFrog.h"

//===================================================================//
//=======================EXPERIMENTAL=CASES===========================//
//===================================================================//
// ==== DamBreak mDBC-reinitialization === //
//#include "../cases/damBreakLong_mDBC/SPH_simulationSymplectic_reinit.h"
//#include "../cases/damBreakLong_mDBC/SPH_simulationVerlet_reinit.h"

// ========= DamBreak GWBC ======== //
//#include "../cases/damBreakLong_GWBC/SPH_simulation.h" /* zatim stale nefunguje */

//===================================================================//

// --- Still Watter BT --- //
//#include "../cases/stillWater_BT/SPH_simulationSymplectic.h"

// --- RALESPH Still Watter --- //
//#include "../cases/ALESPH_stillWater/MUSCL.h"

// --- DamBreak mDBC --- //
//#include "../cases/damBreak_mDBC/SPH_simulationSymplectic.h"

// --- DamBreak BT --- //
//#include "../cases/damBreak_BT/SPH_simulationSymplectic.h"

// --- RSPH DamBreak --- //
//#include "../cases/RSPH_dambreak/SPH_simulation.h"
//#include "../cases/RSPH_dambreak/SPH_simulationSymplectic.h"

// --- RALESPH DamBreak --- //
//#include "../cases/ALESPH_dambreak/MUSCL.h"
//#include "../cases/ALESPH_dambreak/SPH_simulationINT.h"

// --- inflowSlab mDBC --- //
//#include "../cases/inflowSlab_mDBC/SPH_simulationSymplectic.h"
//tohle zatim nefunguje: #include "../cases/inflowSlab_mDBC/SPH_simulationVerlet.h"

//#include "../cases/damBreakLong_BT/SPH_simulationSymplectic.h"
//#include "../cases/damBreakLong_BToprava/SPH_simulationSymplectic.h"
//#include "../cases/ALESPH_dambreakLong/MUSCL.h"

//ALESOON
//===================================================================//
//===========================IN=PROGRESS=============================//
//===================================================================//
//#include "../cases/ALESPH_dambreakLongNEW/SPH_simulationINT.h"
//#include "../cases/ALESPH_dambreakLongNEW/SPH_simulationINTVerlet.h"
//#include "../cases/RSPH_damBreakLong/SPH_simulationSymplectic.h"
//===================================================================//


int main(){

	printf("working... \n");

	SPH_simulation SIMULATION;

	SIMULATION.PREP_SIMULATION_DATA(); //move to simulation
	SIMULATION.INIT();
	SIMULATION.RUN();

	return EXIT_SUCCESS;
}
