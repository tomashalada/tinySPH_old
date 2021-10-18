#pragma once

#include "SPH_defs.h"

/* Artificial viscosity */
real Artificial_Viscosity
(real h, real drs, real drdv, real rho0, real c0, real alpha)
{

	real AV; //Artificial viscosity, \Pi in theory
	real mu; //-> theory

	real visco = MAX(drdv, alpha);

	if (drdv < 0)
	{

		mu = h*drdv / (drs*drs + eps*h*h);
		AV = - visco * c0 * mu / rho0;

	}
	else
	{

		AV = 0;

	}

	return AV;
}

/* Artificial viscosity */
real Laminar_and_SPS_Viscosity
()
{
	real AV = 0;

	return AV;
}
