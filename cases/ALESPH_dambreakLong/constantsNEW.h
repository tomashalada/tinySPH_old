#include "SPH_defs.h"

// CASE: StillWater_mDBC
// NOTE:

// ========== DOMAIN LIMITS ========== //
const real x_0 = 0.; //min point x
const real x_1 = 2.; //max point x
const real y_0 = 0.; //min point y
const real y_1 = 1.; //max point y

// ========== SPH BASIC CONSTANTS ========== //
const real dp = 0.005; //initial particle distance
//const real dp = 0.01; //initial particle distance
const real kap = 2.1; //kernel domain constant (2 for Wendland kernel!)

const real hh = 1.0*sqrt(2*dp*dp); //smoothing length
//const real hh = 2.*dp; //smoothing length

//const real dt = 0.0001; //time step (initial time step)
const real dt = 0.00004; //time step (initial time step)

const	int nvl = std::ceil(kap*hh/dp); //number of boundary layers for ghost/virtual layer based BC

// ========== SPH CONSTANTS ========== //
const	real cs = 44.3; //speed of sound
const real visco = 0.02; //viscosity (coef. of artificial viscosity for WCSPH, viscosity value otherwise)
const	real delta = 0.1; //delta-SPH constant

// ========== SIMULATION CONTROL ========== //
const real t_max = 2.; //NOW ACTIVE YET

unsigned int step = 0; //Time in steps
const unsigned int step_end = 50000; //Endtime in steps;

const unsigned int save_output_interval = 250; //Interval to save results
const unsigned int save_output_interval_interpolated = 2500; //Interval to save interpolated results

// ========== SIMULATION CONTROL ========== //

//Case directory
std::string casePATH = "/home/tomas/Documents/__sovler/tinySPH_double_mr/cases/ALESPH_dambreakLong/";
//std::string casePATH = "/storage/praha1/home/haladto1/tinySPH/damBreakLong_mDBC/";

std::string fileName_resultsFluidOnly = casePATH + "OUTPUT/resultsFluidOnly/partFluidOnly_"; //particles: fluid only
std::string fileName_resultsAll = casePATH + "OUTPUT/resultsAll/partAll_"; //particles: all
std::string fileName_particlesOut = casePATH + "OUTPUT/resultsAll/partAll_"; //particles: out of domain

std::string fileName_resultsInterpol = casePATH + "OUTPUT/resultsInterpolated/fieldInterpolated_"; //interpolated variables

std::string fileName_initCond = casePATH + "OUTPUT/initialCondition.vtk"; //initial condition
std::string fileName_info = casePATH + "OUTPUT/simulationInfo.txt"; //simulation info
std::string fileName_grid = casePATH + "OUTPUT/grid.vtk"; //find neighbour grid (case grid)

// ========== MEASUREMENT CONTROL ========== //
// ---------- water level ------------ //
const unsigned int save_WL_interval = 50; //Interval to save results
realvec wl_h1_pos = {0.3, 0.}; // position of water level h1 sensore
std::string fileName_waterLevel_h1 = casePATH + "OUTPUT/waterLevel_h1.dat"; //water level on h1 position

realvec wl_h1_pos_inv = {0.3, 0.45}; // position of water level h1 sensore
std::string fileName_waterLevel_h1_inv = casePATH + "OUTPUT/waterLevel_h1_inv.dat"; //water level on h1 position

realvec wl_h2_pos = {0.865, 0.}; // position of water level h1 sensore
std::string fileName_waterLevel_h2 = casePATH + "OUTPUT/waterLevel_h2.dat"; //water level on h1 position

realvec wl_h2_pos_inv = {0.865, 0.45}; // position of water level h1 sensore
std::string fileName_waterLevel_h2_inv = casePATH + "OUTPUT/waterLevel_h2_inv.dat"; //water level on h1 position

realvec wl_h3_pos = {1.114, 0.}; // position of water level h1 sensore
std::string fileName_waterLevel_h3 = casePATH + "OUTPUT/waterLevel_h3.dat"; //water level on h1 position

realvec wl_h3_pos_inv = {1.114, 0.45}; // position of water level h1 sensore
std::string fileName_waterLevel_h3_inv = casePATH + "OUTPUT/waterLevel_h3_inv.dat"; //water level on h1 position

realvec wl_h4_pos = {1.3625, 0.}; // position of water level h1 sensore
std::string fileName_waterLevel_h4 = casePATH + "OUTPUT/waterLevel_h4.dat"; //water level on h1 position

realvec wl_h4_pos_inv = {1.3625, 0.45}; // position of water level h1 sensore
std::string fileName_waterLevel_h4_inv = casePATH + "OUTPUT/waterLevel_h4_inv.dat"; //water level on h1 position

// ---------- pressure ------------ //
const unsigned int save_pressure_interval = 15; //Interval to save results
realvec p_h1_pos = {1.61, 0.003}; // position of water level h1 sensore
std::string fileName_waterLevel_p1 = casePATH + "OUTPUT/pressure_p1.dat"; //water level on h1 position

realvec p_h2_pos = {1.61, 0.015}; // position of water level h1 sensore
std::string fileName_waterLevel_p2 = casePATH + "OUTPUT/pressure_p2.dat"; //water level on h1 position

realvec p_h3_pos = {1.61, 0.03}; // position of water level h1 sensore
std::string fileName_waterLevel_p3 = casePATH + "OUTPUT/pressure_p3.dat"; //water level on h1 position

realvec p_h4_pos = {1.61, 0.08}; // position of water level h1 sensore
std::string fileName_waterLevel_p4 = casePATH + "OUTPUT/pressure_p4.dat"; //water level on h1 position

