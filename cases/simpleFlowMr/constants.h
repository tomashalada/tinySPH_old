#include "SPH_defs.h"

	/*** Simulation parameters & settings ***/

	double x_0 = 0.; // temp
	double x_1 = 2.; // temp
	double y_0 = 0.; // temp
	double y_1 = 1.; // temp

	real dpLR = 0.005; //temp
	real dpHR = 0.003; //temp
	//real dpHR = 0.005; //temp

	// real dpLR = 0.005; //temp
	// real dpHR = 0.0025; //temp
	// //real dpHR = 0.005; //temp

	real dt = 0.0001; //temp
	real kap = 1.;

	unsigned int save_output_interval = 25; //temp, default 25
	unsigned int save_output_interval_interpolated = 1000;
	unsigned int step = 0; //Time in steps
	unsigned int step_end = 5001; //Endtime in steps;

 std::string output_file_nameLR = "resultsLR/particles_";
 std::string output_file_nameHR = "resultsHR/particles_";
	//std::string output_file_name_p_out = "p_out/particles_";

 std::string output_file_nameLRfo = "resultsLRfo/particles_";
 std::string output_file_nameHRfo = "resultsHRfo/particles_";

	/*** Inlet and outlet zone settings ***/

	//Inlet constant buffer in LR region.
	realvec LRinlet_vel = {2. , 0.};
	real LRinlet_x0 = 0.1 + 4*dpLR;
	real LRinlet_y0 = 0.0 + dpLR;

	//Outlet dynamic buffer in LR region, with dynamic height
	real LRoutletDyn_x0 = 0.9-4*dpLR;
	real LRoutletDyn_y0 = 0.0 + dpLR;
	real LRoutletDyn_xm = 0.9;
	real LRoutletDyn_ym = 0.2;

	/*** Help & temp variables ***/

	double hhLR = 1*sqrt(2*dpLR*dpLR);
	double hhHR = 1*sqrt(2*dpLR*dpHR);
	//double hh = 0.02;
	int nvlLR = std::ceil(kap*hhLR/dpLR);
	int nvlHR = std::ceil(kap*hhHR/dpHR);

// OUTPUT FILES AND DIRECTORIES
std::string casePATH = "/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/simpleFlowMr/";

std::string fileName_resultsFluidOnlyHR = casePATH + "OUTPUT/resultsFluidOnly/partFluidOnlyHR_";
std::string fileName_resultsFluidOnlyLR = casePATH + "OUTPUT/resultsFluidOnly/partFluidOnlyLR_";

std::string fileName_resultsAllHR = casePATH + "OUTPUT/resultsAll/partAllHR_";
std::string fileName_resultsAllLR = casePATH + "OUTPUT/resultsAll/partAllLR_";

std::string fileName_resultsInterpolHR = casePATH + "OUTPUT/resultsInterpolated/fieldInterpolatedHR_";
std::string fileName_resultsInterpolLR = casePATH + "OUTPUT/resultsInterpolated/fieldInterpolatedLR_";

std::string fileName_initCondHR = casePATH + "OUTPUT/initialConditionHR.vtk";
std::string fileName_initCondLR = casePATH + "OUTPUT/initialConditionLR.vtk";

std::string fileName_infoHR = casePATH + "OUTPUT/simulationInfoHR.txt";
std::string fileName_infoLR = casePATH + "OUTPUT/simulationInfoLR.txt";



