#include "SPH_defs.h"

	/*** Simulation parameters & settings ***/

	double x_0 = 0.; // temp
	double x_1 = 2.; // temp
	double y_0 = 0.; // temp
	double y_1 = 1.; // temp

	real dp = 0.005; //temp
	//real dpHR = 0.005; //temp

	// real dp = 0.005; //temp
	// real dpHR = 0.0025; //temp
	// //real dpHR = 0.005; //temp

	//real dt = 0.0001; //temp
	real dt = 0.00005; //temp
	real kap = 1.;

	unsigned int save_output_interval = 25; //temp, default 25
	unsigned int save_output_interval_interpolated = 1000;
	unsigned int step = 0; //Time in steps
	unsigned int step_end = 10001; //Endtime in steps;

 std::string output_file_name = "results/particles_";
	//std::string output_file_name_p_out = "p_out/particles_";

 std::string output_file_namefo = "resultsfo/particles_";

	/*** Inlet and outlet zone settings ***/

	//Inlet constant buffer in  region.
	realvec inlet_vel = {2. , 0.};
	real inlet_x0 = 0.1 + 4*dp;
	real inlet_y0 = 0.0 + dp;

	//Outlet dynamic buffer in  region, with dynamic height
	real outletDyn_x0 = 0.9-4*dp;
	real outletDyn_y0 = 0.0 + dp;
	real outletDyn_xm = 0.9;
	real outletDyn_ym = 0.2;

	/*** Help & temp variables ***/

	double hh = 1.1*sqrt(2*dp*dp);
	//double hh = 0.02;
	int nvl = std::ceil(kap*hh/dp);

// OUTPUT FILES AND DIRECTORIES
std::string casePATH = "/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/simpleFlowSr_mdbcLattice/";

std::string fileName_resultsFluidOnly = casePATH + "OUTPUT/resultsFluidOnly/partFluidOnly_";

std::string fileName_resultsAll = casePATH + "OUTPUT/resultsAll/partAll_";

std::string fileName_resultsInterpol = casePATH + "OUTPUT/resultsInterpolated/fieldInterpolatedTESTsFINE_";

std::string fileName_initCond = casePATH + "OUTPUT/initialCondition.vtk";

std::string fileName_info = casePATH + "OUTPUT/simulationInfo.txt";



