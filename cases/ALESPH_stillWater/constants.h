#include "SPH_defs.h"

// CASE: StillWater_mDBC
// NOTE:

// ========== DOMAIN LIMITS ========== //
const real x_0 = 0.; //min point x
const real x_1 = 1.; //max point x
const real y_0 = 0.; //min point y
const real y_1 = 1.; //max point y

// ========== SPH BASIC CONSTANTS ========== //
const real dp = 0.02; //initial particle distance
const real kap = 2.1; //kernel domain constant (2 for Wendland kernel!)

//const real hh = 1.0*sqrt(2*dp*dp); //smoothing length
const real hh = 2.*dp; //smoothing length

const real dt = 0.0002; //time step (initial time step)

const	int nvl = std::ceil(kap*hh/dp); //number of boundary layers for ghost/virtual layer based BC

// ========== SPH CONSTANTS ========== //
const	real cs = 42.48; //speed of sound
const real visco = 0.01; //viscosity (coef. of artificial viscosity for WCSPH, viscosity value otherwise)
const	real delta = 0.1; //delta-SPH constant


// ========== SIMULATION CONTROL ========== //
const real t_max = 4.; //NOW ACTIVE YET

unsigned int step = 0; //Time in steps
const unsigned int step_end = 20001; //Endtime in steps;

const unsigned int save_output_interval = 200; //Interval to save results
const unsigned int save_output_interval_interpolated = 1000; //Interval to save interpolated results

// ========== SIMULATION CONTROL ========== //

//Case directory
std::string casePATH = "/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/ALESPH_stillWater/";

std::string fileName_resultsFluidOnly = casePATH + "OUTPUT/resultsFluidOnly/partFluidOnly_"; //particles: fluid only
std::string fileName_resultsAll = casePATH + "OUTPUT/resultsAll/partAll_"; //particles: all
std::string fileName_particlesOut = casePATH + "OUTPUT/resultsAll/partAll_"; //particles: out of domain

std::string fileName_resultsInterpol = casePATH + "OUTPUT/resultsInterpolated/fieldInterpolated_"; //interpolated variables

std::string fileName_initCond = casePATH + "OUTPUT/initialCondition.vtk"; //initial condition
std::string fileName_info = casePATH + "OUTPUT/simulationInfo.txt"; //simulation info
std::string fileName_grid = casePATH + "OUTPUT/grid.vtk"; //find neighbour grid (case grid)



