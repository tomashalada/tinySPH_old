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

real densRiemannLinearizedWithLimiter
(real rhoL, real rhoR, real vL, real vR, real rhoAvg, real c)
{

	real rho;

	real p;
	real eta = 3.;
	real beta = MIN(eta*MAX(vL-vR, 0), c);
	rho = 0.5*(rhoL + rhoR) + 0.5*beta*(vL - vR)*rhoAvg/(c*c);
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

real pRiemannLinearizedNecoNeco
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

//help function
real cs
(real rho, real rho0, real c0)
{

	real cs;

	real gamma = 2.;
	cs = c0*pow((rho/rho0), 0.5*(gamma-1));

	return cs;

}

real minmod
(real r)
{

	real fi;
	fi = MAX(0, MIN(1,r));
	return fi;

}

//MUSCL: Update L and R state
real MUSCLleft
(real A, realvec gradA, realvec dr)
{

	real AL = A + 0.5*(gradA.x * dr.x + gradA.y * dr.y);

	return AL;

}

real MUSCLright
(real A, realvec gradA, realvec dr)
{

	real AR = A - 0.5*(gradA.x * dr.x + gradA.y * dr.y);

	return AR;

}


//MUSCL: Update L and R state with limiter
real MUSCLleft
(real A, realvec gradA, realvec dr, real alpha)
{

	real AL = A - alpha*0.5*(gradA.x * dr.x + gradA.y * dr.y);

	return AL;

}

real MUSCLright
(real A, realvec gradA, realvec dr, real alpha)
{

	real AR = A + alpha*0.5*(gradA.x * dr.x + gradA.y * dr.y);

	return AR;

}


