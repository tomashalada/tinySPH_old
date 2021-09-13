#include "SPH_defs.h"


// struct SIMULATION_PARAMETERS
// {
//
// 	/*** Main, global parameters. ***/
// 	real dpLR; //initial particle spacing LR region
// 	real dpHR; //initial particle spacing HR region
// 	real dt; //time step
// 	real kap; //kernel cut-off distance constant
//
//
// 	/*** Time and output settings. ***/
//  double t;
//  double t_fin;
// 	unsigned int step_end; //endime in steps
// 	unsigned int save_output_interval; //output interval for the results
//
// 	/*** Output files settings. ***/
// 	std::string output_file_nameLR; //move
// 	std::string output_file_nameHR; //move
// 	//std::string output_file_name_p_out = "p_out/particles_"; //move
//
// 	/*** Inlet and outlet parameters. ***/
// 	realvec inlet_vel;
// 	realvec outlet_vel;
// 	unsigned int outlet_type; //0 for fixed velocity, 1 for extrapolated velocity
//
//
// 	/* inlet position */
//
// 	/* outlet position */
//
// };


/*** Main, global parameters. ***/
extern real _dpLR = 0.01; //initial particle spacing LR region
extern real _dpHR = 0.0075; //initial particle spacing HR region
extern real _dt = 0.0001; //time step
extern real _kap = 1.; //kernel cut-off distance constant


/*** Time and output settings. ***/
extern double t;
extern double t_fin;
extern unsigned int step_end = 5001; //endime in steps
extern unsigned int save_output_interval = 25; //output interval for the results

/*** Output files settings. ***/
std::string output_file_nameLR = "resultsLR/particles_"; //move
std::string output_file_nameHR = "resultsHR/particles_"; //move
//std::string output_file_name_p_out = "p_out/particles_"; //move

/*** Inlet and outlet parameters. ***/
realvec inlet_vel = {2., 0.};
realvec outlet_vel = {0., 0.};
unsigned int outlet_type = 0; //0 for fixed velocity, 1 for extrapolated velocity


/* inlet position */

/* outlet position */




