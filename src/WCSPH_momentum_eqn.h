#include "SPH_defs.h"

real p_temp = 0; //pressure term
real visco = 0; //viscosity term

//Assign acceleration terms
//p_temp = ap/pow(arho, 2) + np/pow(nrho,2);
p_temp = (ap + np)/(arho*nrho);
visco = Artificial_Viscosity(h, drs, drdv, 0.5*(arho + nrho), c0, particles.data_const.avisc);
