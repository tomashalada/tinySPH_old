#pragma once

#include "SPH_defs.h"
#include "SPH_pressure.h"

// ========== ARTIFICIAL VISCOSITY ========== //
real Artificial_Viscosity
(real h, real drs, real drdv, real rho0, real c0, real alpha)
{

	real AV = 0; //Artificial viscosity, \Pi in theory
	real mu; //-> theory

	real visco = MAX(drdv, alpha);

	if (drdv < 0)
	{

		mu = h*drdv / (drs*drs + eps*h*h);
		AV = - visco * c0 * mu / rho0;

	}

	return AV;
}

// ========== DENS. DIFFUSION TERMS ========== //
real Density_Diffusion_term_MOLTENI
(real h, real drs, real drdW, real arho, real nrho, real c0, real delta, real m)
{

	real DT;
	real PSI_s; //psi scalar (without dr)

	PSI_s = 2. * c0 * h * delta * (nrho - arho) / (drs*drs + eps*h*h);
	DT = PSI_s * drdW * m / nrho;

	//PSI_s = 2. * c0 * h * delta * (arho/nrho - 1.) / (drs*drs + eps*h*h);
	//DT = PSI_s * drdW * m ;

	return DT;

}

real Density_Diffusion_term_FOURTAKAS
(real h, real drs, real dr_z, real drdW, real arho, real nrho, real rho0, real c0, real delta, real m, real gravity)
{

	real DT;
	real PSI_s; //psi scalar (without dr)

	const real cb = c0 * c0 * rho0 / 7.;
	const real dpH = 1. + dr_z * rho0 * gravity / cb;
	const real drhoH = rho0 * pow(dpH, 1./7) - rho0;

	PSI_s = 2. * c0 * h * delta * ((nrho - arho) - drhoH) / (drs*drs + eps*h*h);
	DT = PSI_s * drdW * m / nrho;

	return DT;

}

