#pragma once
#include "SPH_defs.h"
#include "SPH_pressure.h"

std::array<real, 3> F
(real rho, real u, real v, real p)
{

	std::array<real, 3> F;

	F[0] = rho*u + rho*v;
	F[1] = rho*u*u + p + rho*u*v;
	F[2] = rho*u*v + rho*v*v +p;

	return F;
}

std::array<real, 3> F
(std::array<real,3> W, real rho0, real cs)
{

	std::array<real, 3> F;
	const real p = Compute_Pressure(W[0], rho0, cs);

	F[0] = W[1] + W[2];
	F[1] = W[1]*W[1]/W[0] + p + W[1]*W[2]/W[0];
	F[2] = W[1]*W[2]/W[0] + W[2]*W[2]/W[0] + p;

	return F;
}


real MUSCL_superBee_L
(real A, realvec gradF, realvec gradB, realvec dr)
{
	const real beta = 2.;

	real AL; //left state

	real D;
	real DF = gradF.x * dr.x  + gradF.y * dr.y;
	real DB = gradB.x * dr.x  + gradB.y * dr.y;

	if(DF > 0)
	{ D = MAX(0, MAX(MIN(beta*DB,DF), MIN(DB,DF))); }
	else if (DF < 0)
	{ D = MIN(0, MIN(MAX(beta*DB,DF), MAX(DB,DF))); }
	else
	{ D = 0.; }

	AL = A - D;

	return AL;
}

real MUSCL_superBee_R
(real A, realvec gradF, realvec gradB, realvec dr)
{
	const real beta = 2.;

	real AL; //left state

	real D;
	real DF = gradF.x * dr.x  + gradF.y * dr.y;
	real DB = gradB.x * dr.x  + gradB.y * dr.y;

	if(DF > 0)
	{ D = MAX(0, MAX(MIN(beta*DB,DF), MIN(DB,DF))); }
	else if (DF < 0)
	{ D = MIN(0, MIN(MAX(beta*DB,DF), MAX(DB,DF))); }
	else
	{ D = 0.; }

	AL = A + D;

	return AL;
}
std::array<real, 3> MUSCL_final
(std::array<real, 3> WaL, std::array<real, 3> WaR, std::array<real, 3> aW, real dt, real dx, real rho0, real cs)
{
	std::array<real, 3> W;

	W[0] = aW[0] + 0.5*dt/dx*(F(WaL,rho0,cs)[0] - F(WaL,rho0,cs)[0]);
	W[1] = aW[1] + 0.5*dt/dx*(F(WaL,rho0,cs)[1] - F(WaL,rho0,cs)[1]);
	W[2] = aW[2] + 0.5*dt/dx*(F(WaL,rho0,cs)[2] - F(WaL,rho0,cs)[2]);

	return W;
}
