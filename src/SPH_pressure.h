#pragma once

#include "SPH_defs.h"

real Compute_Pressure
(real rho, real rho0, real c0)
{

	real p;
	//c0 = 40; //temp, tst

	real gamma = 7.; //Weakly compressible state eq. Poisson constant
	real b_const; //Weakly compressible state eq. constant
	real rho_r; //Weakly compressible state eq. relative density

	b_const = c0*c0*rho0/gamma;
	rho_r = rho/rho0;

	/* Debug */
	//std::cout << "PRESSURE -> b0 constant: " << b_const << std::endl;

	/* Weakly compressible state equation. */
	p = b_const*( pow(rho_r, gamma) - 1 );



	return p;

}


real Compute_Pressure2
(real rho, real rho0, real c0)
{

	real p;
	//c0 = 40; //temp, tst

	real gamma = 7.; //Weakly compressible state eq. Poisson constant
	real b_const; //Weakly compressible state eq. constant
	real rho_r; //Weakly compressible state eq. relative density

	b_const = c0*c0*rho0/gamma;
	rho_r = rho/rho0;

	/* Debug */
	//std::cout << "PRESSURE -> b0 constant: " << b_const << std::endl;

	/* Weakly compressible state equation. */
	//p = b_const*( pow(rho_r, gamma) - 1 );
	p =c0*c0*( rho - rho0 );
	//p =c0*c0*( rho0 - rho );



	return p;

}
