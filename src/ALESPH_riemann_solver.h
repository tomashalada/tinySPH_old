#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

real densRiemannLinearized
(real rhoL, real rhoR, real vL, real vR, real rhoAvg, real c)
{

	real rho;
	rho = 0.5*(rhoL + rhoR) + 0.5*(vL - vR)*rhoAvg/c;
	//rho = 0.5*(rhoL + rhoR);// - 0.5*(vR - vL)*rhoAvg/c;

	return rho;

}

real pRiemannLinearized
(real rhoL, real rhoR, real vL, real vR, real pL, real pR, real rhoAvg, real c)
{

	real p;
	//p = 0.5*(pL + pR) - 0.5*(vR - vL)*rhoAvg*c;
	//p = 0.5*(pL + pR) + 0.5*(vL - vR)*rhoAvg*c;
	p = 0.5*(pL + pR) + 0.5*(vL - vR)*rhoAvg*c;

	/// p = 0.5*(pL + pR) ;//+ 0.5*(vL - vR)*rhoAvg*c;

	return p;

}

real pRiemannLinearizedWithLimiter
(real rhoL, real rhoR, real vL, real vR, real pL, real pR, real rhoAvg, real c)
{

	real p;
	real eta = 3.;
	real beta = MIN(eta*MAX(vL-vR, 0), c);

	//p = 0.5*(pL + pR) - 0.5*(vR - vL)*rhoAvg*c;
	//p = 0.5*(pL + pR) + 0.5*(vL - vR)*rhoAvg*c;

	///FUNGUJE: p = 0.5*(pL + pR) + 0.5*(vL - vR)*rhoAvg*c;
	p = 0.5*(pL + pR) + 0.5*(vL - vR)*rhoAvg*beta;

	/// p = 0.5*(pL + pR) ;//+ 0.5*(vL - vR)*rhoAvg*c;

	return p;

}

real velRiemannLinearized
(real rhoL, real rhoR, real vL, real vR, real rhoAvg, real c)
{

	real v;
	v = 0.5*(vR + vL) - 0.5*(rhoR - rhoL)*c/rhoAvg;
	//v = 0.5*(vR + vL);// - 0.5*(rhoR - rhoL)*c/rhoAvg;

	return v;

}

real velRiemannLinearizedwithPressure
(real rhoL, real rhoR, real vL, real vR, real pL, real pR, real rhoAvg, real c)
{

	real v;
	//v = 0.5*(vR + vL) - 0.5*(pL - pR)*c*rhoAvg;
	//v = 0.5*(vR + vL) - 0.5*(pL - pR)*c*rhoAvg;
	v = 0.5*(vR + vL) + 0.5*(pL - pR)/(c*rhoAvg);

	/// v = 0.5*(vR + vL);// + 0.5*(pL - pR)/(0.5*c*rhoAvg);

	return v;

}
