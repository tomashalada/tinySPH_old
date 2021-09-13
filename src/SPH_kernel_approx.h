#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

/* Approximate velocity with kernel aproximation to certain point */
realvec Kernel_velocity_approximation
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	realvec VEL_APROX = {0. , 0.};

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;

	int npx, npy;
	int ncx, ncy;

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	realvec nv;
	realvec dr;
	real drs;
	real nrho;
	//real drdv;

	double *kernel;

	realvec vWsum = {0., 0.};
	real Wsum = 0;
	realvec vvv = {0., 0.};
	//realvec vWsum;
	//real Wsum ;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		vWsum = vWsum + nv*kernel[0]; // += not overloaded yet
		vvv = vvv + nv*kernel[0]*(m/nrho);

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//VEL_APROX = VEL_APROX + vWsum; //??? make it better
	//VEL_APROX = VEL_APROX + vvv; //??? make it better

	if(Wsum == 0){VEL_APROX = {0. ,0.};}

	return VEL_APROX;
}

/* Approximate velocity with kernel aproximation to certain point */
real Kernel_density_approximation
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	real dens = 1000.;

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	real nrho;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	real rhoWsum = 0;
	real Wsum = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += nrho*kernel[0]; // += not overloaded yet


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	dens =  rhoWsum / Wsum; //??? make it better
	if(Wsum == 0){dens = 1000;}

	return dens;
}


/* Approximate velocity with kernel aproximation to certain point */
real Kernel_get_waterlevel
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	//real dens = 1000.;
	real am = 0;

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;
	real dp = particles.data_const.dp;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr = {0., 0};
	real nrho = 0;
	realvec dr = {0., 0};
	real drs = 0;
	//real drdv;

	double *kernel;

	real mWsum = 0;
	real Wsum = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		/* I think mass sumation should be better bacause some strange stuff near outlet buffer in density field*/
		//mWsum += dp*dp*nrho*kernel[0]; // += not overloaded yet
		mWsum += m*kernel[0]; // += not overloaded yet


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	//am =  mWsum / Wsum; //??? make it better
	am =  mWsum; //??? make it better

	return am;
}
/*

real get_reference_mass_value
(Particle_system &particles)
{

	real mfull;

	return mfull;

}
*/

real Water_Elevation
(Particle_system &particles, Simulation_data simulation_data, real x, real y0, real dp)//position, where we want to know water elevation
{

 real ELEVATION = 0;
	/* get reference mass value */
	real m = 0;
	real am = 0;
	realvec pos = {0., 0.};
	pos.x = x;

	realvec posref = {0.2 , 0.1};
	//real mfull = get_reference_mass_value(particles);
	real mfull = Kernel_get_waterlevel(particles, simulation_data, posref);
	//real mfull = particles.data_const.m;

	/* Debug */
	std::cout << "WATER-ELEVATION: Full mass: " << mfull << " reference coordinates: [" << posref.x << "," << posref.y << "]" << std::endl;
	//mfull = Kernel_get_waterlevel(particles, simulation_data, posref);
	mfull = 1000.;


	for(int i = 0; i < int(0.5/dp); i++){

		pos.y = y0 + i*(dp/4); //haha
	 am = Kernel_get_waterlevel(particles,simulation_data, pos);

		/* Debug */
		//std::cout << "WATER-ELEVATION: Actuall mass: " << am << " coordinates: [" << pos.x << "," << pos.y << "]" << std::endl;

		if(am < mfull/2 || std::isnan(am))
		{
		ELEVATION = pos.y;
		//std::cout << "WATER-ELEVATION: Actuall mass: " << am << " coordinates: [" << pos.x << "," << pos.y << "]" << std::endl;
		break;
		}


	}

		/* Debug */
		std::cout << "WATER-ELEVATION: ELEVATION " << am << " coordinates: [" << pos.x << "," << pos.y << "] is " << ELEVATION << std::endl;
	return ELEVATION;

}

//------------------------------------------------------------------------------------------
/* Approximate velocity with kernel aproximation to certain point */
realvec Kernel_velocity_approximationTRZ
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	realvec VEL_APROX = {0. , 0.};

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	realvec nv;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	realvec vWsum = {0., 0.};
	real Wsum = 0;
	//realvec vWsum;
	//real Wsum ;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];

		/* Debug */
		//std::cout << "KERNEL APROX. TRANSITION ZONE: Active. nv: [" << nv.x <<","<< nv.y << "] on position: [" << nr.x <<","<< nr.y << "]" << std::endl;

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		vWsum = vWsum + nv*kernel[0]; // += not overloaded yet

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//VEL_APROX = VEL_APROX + vWsum; //??? make it better
	if(Wsum == 0){VEL_APROX = {0. ,0.};}

	return VEL_APROX;
}

/* Approximate velocity with kernel aproximation to certain point */
real Kernel_density_approximationTRZ
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	real dens = 1000.;

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	real nrho;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	real rhoWsum = 0;
	real Wsum = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += nrho*kernel[0]; // += not overloaded yet


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	dens =  rhoWsum / Wsum; //??? make it better
	if(Wsum == 0){dens = 1000;}

	return dens;
}

/* Approximate velocity with kernel aproximation to certain point */
real Kernel_density_approximation_mDBC
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	real dens = 1000.;

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	real nrho;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	real rhoWsum = 0;
	real Wsum = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += nrho*kernel[0]; // += not overloaded yet


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	dens =  rhoWsum / Wsum; //??? make it better
	if(Wsum == 0){dens = 1000;}

	return dens;
}

//*************************************************************************************************//
real determinant(real a11, real a12, real a13,
																	real a21, real a22, real a23,
																	real a31, real a32, real a33)
{

	//return a11*a22*a33 + a21*a32*a13 + a31*a21*a23 - a13*a22*a31 - a23*a32*a11 - a33*a12*a21;
	//std::cout << "----> DET DIF: " << (a11*a22*a33 + a21*a32*a13 + a31*a12*a23 - a13*a22*a31 - a23*a32*a11 - a33*a12*a21) - (a11*(a22*a33 - a23*a32) + a21*(a23*a13 - a12*a33) + a31*(a12*a23 - a22*a13)) << std::endl;
	return (a11*(a22*a33 - a23*a32) + a21*(a23*a13 - a12*a33) + a31*(a12*a23 - a22*a13));

}


/* Approximate velocity with kernel aproximation to certain point */
realvec Kernel_velocity_approximation_MATRIX
(Particle_system &particles, Simulation_data simulation_data, realvec ar, realvec gnd)
{

	realvec VEL_APROX = {0. , 0.}; //vel. to return
	realvec VEL_APROX_temp = {0. , 0.};

	/* Corection matrix, velocity sum., velocity gradient sum. */

	real a11x, a11inv;
	real a12x, a12inv;
	real a13x, a13inv;

	real a21x, a21inv;
	real a22x, a22inv;
	real a23x, a23inv;

	real a31x, a31inv;
	real a32x, a32inv;
	real a33x, a33inv;

	real dudx, dudy;
	real dvdx, dvdy;

	a11x = 0, a11inv = 0;
	a12x = 0, a12inv = 0;
	a13x = 0, a13inv = 0;

	a21x = 0, a21inv = 0;
	a22x = 0, a22inv = 0;
	a23x = 0, a23inv = 0;

	a31x = 0, a31inv = 0;
	a32x = 0, a32inv = 0;
	a33x = 0, a33inv = 0;

	realvec vs = {0., 0.};

	dudx = 0, dudy = 0;
	dvdx = 0, dvdy = 0;


	/* Get SPH constants. */

	real kh, of, h, m;

	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	h = particles.data_const.h;
	m = particles.data_const.m;

	/* Get cell and nb's cells for point of aprotixmation. */

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	/* Interaction data */

	realvec nr = {0., 0.};
	realvec nv = {0., 0.};
	realvec dr = {0., 0.};
	real drs = 0;
	real nrho = 0;

	double *kernel;

	realvec vWsum = {0., 0.};
	real Wsum = 0;

	real vWx = 0;
	real vWy = 0;


	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		//vWsum = vWsum + nv*kernel[0]* (m/nrho); // += not overloaded yet

		vWx += nv.x * kernel[0] * (m/nrho);
		vWy += nv.y * kernel[0] * (m/nrho);


		dudx += nv.x*dr.x*kernel[1] * (m/nrho);
		dudy += nv.x*dr.y*kernel[1] * (m/nrho);
		dvdx += nv.y*dr.x*kernel[1] * (m/nrho);
		dvdy += nv.y*dr.y*kernel[1] * (m/nrho);

		a11x += kernel[0]*(m/nrho);
		a12x += dr.x*kernel[0]*(m/nrho);
		a13x += dr.y*kernel[0]*(m/nrho);

		a21x += dr.x*kernel[1]*(m/nrho);
		a22x += dr.x*dr.x*kernel[1]*(m/nrho);
		a23x += dr.x*dr.y*kernel[1]*(m/nrho);

		a31x += dr.y*kernel[1]*(m/nrho);
		a32x += dr.y*dr.x*kernel[1]*(m/nrho);
		a33x += dr.y*dr.y*kernel[1]*(m/nrho);

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;

		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	real det;
	det = determinant(a11x, a12x, a13x, a21x, a22x, a23x, a31x, a32x, a33x);

	//std::cout << "[ KERNEL GHOST NODE AROXIMATION ] Determinant: " << det << " vW = [" << vWx << ","<< vWy << "]"<< std::endl;
	if(det > 0.001)
	{
		a11inv =  (a22x*a33x-a23x*a32x)/det;
		a12inv = -(a12x*a33x-a13x*a32x)/det;
		a13inv =  (a12x*a23x-a13x*a22x)/det;

		a21inv = -(a21x*a33x-a23x*a31x)/det;
		a22inv =  (a11x*a33x-a13x*a31x)/det;
		a23inv = -(a11x*a23x-a13x*a21x)/det;

		a31inv =  (a21x*a32x-a22x*a31x)/det;
		a32inv = -(a11x*a32x-a12x*a31x)/det;
		a33inv =  (a11x*a22x-a12x*a21x)/det;

		/*
		std::cout << "[KERNEL APROX.] Matrix | Inv(Matrix):" << std::endl;
		std::cout << a11x << " " << a12x << " " << a13x << "  | " << a11inv << " " << a12inv << " " << a13inv << std::endl;
		std::cout << a21x << " " << a22x << " " << a23x << "  | " << a21inv << " " << a22inv << " " << a23inv << std::endl;
		std::cout << a31x << " " << a32x << " " << a33x << "  | " << a31inv << " " << a32inv << " " << a33inv << std::endl;
		std::cout << "[KERNEL APROX.] Matrix [" << a11x << ", " << a12x << ", " << a13x << "; " \
																																										<< a21x << ", " << a22x << ", " << a23x << "; " \
																																										<< a31x << ", " << a32x << ", " << a33x << "]" << std::endl;

		std::cout << "[KERNEL APROX.] Inv(Matrix) [" << a11inv << ", " << a12inv << ", " << a13inv << "; " \
																																															<< a21inv << ", " << a22inv << ", " << a23inv << "; " \
																																															<< a31inv << ", " << a32inv << ", " << a33inv << "]" << std::endl;

		*/



	//Compute ghost note velocity
	VEL_APROX_temp.x = vWx * a11inv + dudx*a12inv + dudy*a13inv;
	VEL_APROX_temp.y = vWy * a11inv + dvdx*a12inv + dvdy*a13inv;

	real a11 = -(vWx*a21inv + dudx*a22inv + dudy*a23inv);
	real a12 = -(vWy*a21inv + dvdx*a22inv + dvdy*a23inv);
	real a21 = -(vWx*a31inv + dudx*a32inv + dudy*a33inv);
	real a22 = -(vWy*a31inv + dvdx*a32inv + dvdy*a33inv);

	VEL_APROX.x =( VEL_APROX_temp.x + a11*gnd.x + a21*gnd.y);
	VEL_APROX.y =( VEL_APROX_temp.y + a12*gnd.x + a22*gnd.y);

	std::cout << "[KERNEL APROX.] VEL_APROX_temp [" << VEL_APROX_temp.x << ", " << VEL_APROX_temp.y << "]" << std::endl;
	std::cout << "[KERNEL APROX.] FitMatrix [" << a11 << ", " << a12 << ", " << a21 << "; " << a22 <<  "]" << std::endl;
	std::cout << "[KERNEL APROX.] VEL_APROX [" << VEL_APROX.x << ", " << VEL_APROX.y << "]" << std::endl;

	//std::cout << " --> Case A: ACTIVE." << std::endl;

	}
	else if (a11x>0)
	{

	VEL_APROX.x = vWx/a11x;
	VEL_APROX.y = vWy/a11x;

	//// VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//// //VEL_APROX = VEL_APROX + vWsum; //??? make it better

	//std::cout << " --> Case B: ACTIVE." << std::endl;

	}
	else
	{

		VEL_APROX = {0., 0.};

	//std::cout << " --> Case C: ACTIVE." << std::endl;

	}

	return VEL_APROX;


}


/* Approximate velocity with kernel aproximation to certain point */
real Kernel_density_approximation_MATRIX
(Particle_system &particles, Simulation_data simulation_data, realvec ar, realvec gnd)
{

	real dens = 0; //vel. to return
	real dens_temp = 0;

	/* Corection matrix, velocity sum., velocity gradient sum. */

	real a11x, a11inv;
	real a12x, a12inv;
	real a13x, a13inv;

	real a21x, a21inv;
	real a22x, a22inv;
	real a23x, a23inv;

	real a31x, a31inv;
	real a32x, a32inv;
	real a33x, a33inv;

	real drhodx, drhody;

	a11x = 0, a11inv = 0;
	a12x = 0, a12inv = 0;
	a13x = 0, a13inv = 0;

	a21x = 0, a21inv = 0;
	a22x = 0, a22inv = 0;
	a23x = 0, a23inv = 0;

	a31x = 0, a31inv = 0;
	a32x = 0, a32inv = 0;
	a33x = 0, a33inv = 0;

	realvec vs = {0., 0.};

	drhodx = 0, drhody = 0;


	/* Get SPH constants. */

	real kh, of, h, m;

	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	h = particles.data_const.h;
	m = particles.data_const.m;

	/* Get cell and nb's cells for point of aprotixmation. */

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	/* Interaction data */

	realvec nr;
	realvec nv;
	realvec dr;
	real drs;
	real nrho;

	double *kernel;

	real Wsum = 0;
	real rhoWsum = 0;


	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += m*kernel[0];

		drhodx += dr.x*kernel[1] * m;
		drhody += dr.y*kernel[1] * m;

		a11x += kernel[0]*(m/nrho);
		a12x += dr.x*kernel[0]*(m/nrho);
		a13x += dr.y*kernel[0]*(m/nrho);

		a21x += dr.x*kernel[1]*(m/nrho);
		a22x += dr.x*dr.x*kernel[1]*(m/nrho);
		a23x += dr.x*dr.y*kernel[1]*(m/nrho);

		a31x += dr.y*kernel[1]*(m/nrho);
		a32x += dr.y*dr.x*kernel[1]*(m/nrho);
		a33x += dr.y*dr.y*kernel[1]*(m/nrho);

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;

		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	real det;
	det = determinant(a11x, a12x, a13x, a21x, a22x, a23x, a31x, a32x, a33x);

	//std::cout << "[ KERNEL GHOST NODE AROXIMATION ] Determinant: " << det << std::endl;
	//if(det > 0.001)
	if(det > 1000)
	{
		a11inv = (a22x*a33x-a23x*a32x)/det;
		a12inv = -(a12x*a33x-a13x*a32x)/det;
		a13inv = (a12x*a23x-a13x*a22x)/det;

		a21inv = -(a21x*a33x-a23x*a31x)/det;
		a22inv = (a11x*a33x-a13x*a31x)/det;
		a23inv = -(a11x*a23x-a13x*a21x)/det;

		a31inv = (a21x*a32x-a22x*a31x)/det;
		a32inv = -(a11x*a32x-a12x*a31x)/det;
		a33inv = (a11x*a22x-a12x*a21x)/det;

		/*
		std::cout << "[KERNEL APROX.] Matrix | Inv(Matrix):" << std::endl;
		std::cout << a11x << " " << a12x << " " << a13x << "  | " << a11inv << " " << a12inv << " " << a13inv << std::endl;
		std::cout << a21x << " " << a22x << " " << a23x << "  | " << a21inv << " " << a22inv << " " << a23inv << std::endl;
		std::cout << a31x << " " << a32x << " " << a33x << "  | " << a31inv << " " << a32inv << " " << a33inv << std::endl;
		std::cout << "[KERNEL APROX.] Matrix [" << a11x << ", " << a12x << ", " << a13x << "; " \
																																										<< a21x << ", " << a22x << ", " << a23x << "; " \
																																										<< a31x << ", " << a32x << ", " << a33x << "]" << std::endl;

		std::cout << "[KERNEL APROX.] Inv(Matrix) [" << a11inv << ", " << a12inv << ", " << a13inv << "; " \
																																															<< a21inv << ", " << a22inv << ", " << a23inv << "; " \
																																															<< a31inv << ", " << a32inv << ", " << a33inv << "]" << std::endl;

		*/



	//Compute ghost note velocity
	dens_temp = rhoWsum*a11inv + drhodx*a12inv + drhody*a13inv;

	real grx = rhoWsum*a21inv + drhodx*a22inv + drhody*a23inv;
	real gry = rhoWsum*a31inv + drhodx*a32inv + drhody*a33inv;

	dens =( dens_temp + grx*gnd.x + gry*gnd.y);

	//std::cout << " --> Case A: ACTIVE." << std::endl;

	}
	else if (a11x>0)
	{

	dens = rhoWsum/a11x;

	//// VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//// //VEL_APROX = VEL_APROX + vWsum; //??? make it better

	//std::cout << " --> Case B: ACTIVE." << std::endl;

	}
	else
	{

		dens = 1000;

	//std::cout << " --> Case C: ACTIVE." << std::endl;

	}

	return dens;


}

/* Approximate velocity with kernel aproximation to certain point */
real Kernel_density_approximation_MATRIX_mDBC
(Particle_system &particles, Simulation_data simulation_data, realvec ar, realvec gnd)
{

	real dens = 0; //vel. to return
	real dens_temp = 0;

	/* Corection matrix, velocity sum., velocity gradient sum. */

	real a11x, a11inv;
	real a12x, a12inv;
	real a13x, a13inv;

	real a21x, a21inv;
	real a22x, a22inv;
	real a23x, a23inv;

	real a31x, a31inv;
	real a32x, a32inv;
	real a33x, a33inv;

	real drhodx, drhody;

	a11x = 0, a11inv = 0;
	a12x = 0, a12inv = 0;
	a13x = 0, a13inv = 0;

	a21x = 0, a21inv = 0;
	a22x = 0, a22inv = 0;
	a23x = 0, a23inv = 0;

	a31x = 0, a31inv = 0;
	a32x = 0, a32inv = 0;
	a33x = 0, a33inv = 0;

	realvec vs = {0., 0.};

	drhodx = 0, drhody = 0;


	/* Get SPH constants. */

	real kh, of, h, m;

	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	h = particles.data_const.h;
	m = particles.data_const.m;

	/* Get cell and nb's cells for point of aprotixmation. */

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	/* Interaction data */

	realvec nr;
	realvec nv;
	realvec dr;
	real drs;
	real nrho;

	double *kernel;

	real Wsum = 0;
	real rhoWsum = 0;


	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += m*kernel[0];

		drhodx += dr.x*kernel[1] * m;
		drhody += dr.y*kernel[1] * m;

		a11x += kernel[0]*(m/nrho);
		a12x += dr.x*kernel[0]*(m/nrho);
		a13x += dr.y*kernel[0]*(m/nrho);

		a21x += dr.x*kernel[1]*(m/nrho);
		a22x += dr.x*dr.x*kernel[1]*(m/nrho);
		a23x += dr.x*dr.y*kernel[1]*(m/nrho);

		a31x += dr.y*kernel[1]*(m/nrho);
		a32x += dr.y*dr.x*kernel[1]*(m/nrho);
		a33x += dr.y*dr.y*kernel[1]*(m/nrho);

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;

		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	real det;
	det = determinant(a11x, a12x, a13x, a21x, a22x, a23x, a31x, a32x, a33x);

	//std::cout << "[ KERNEL GHOST NODE AROXIMATION ] Determinant: " << det << std::endl;
	if(det > 0.001)
	{
		a11inv = (a22x*a33x-a23x*a32x)/det;
		a12inv = -(a12x*a33x-a13x*a32x)/det;
		a13inv = (a12x*a23x-a13x*a22x)/det;

		a21inv = -(a21x*a33x-a23x*a31x)/det;
		a22inv = (a11x*a33x-a13x*a31x)/det;
		a23inv = -(a11x*a23x-a13x*a21x)/det;

		a31inv = (a21x*a32x-a22x*a31x)/det;
		a32inv = -(a11x*a32x-a12x*a31x)/det;
		a33inv = (a11x*a22x-a12x*a21x)/det;

		/*
		std::cout << "[KERNEL APROX.] Matrix | Inv(Matrix):" << std::endl;
		std::cout << a11x << " " << a12x << " " << a13x << "  | " << a11inv << " " << a12inv << " " << a13inv << std::endl;
		std::cout << a21x << " " << a22x << " " << a23x << "  | " << a21inv << " " << a22inv << " " << a23inv << std::endl;
		std::cout << a31x << " " << a32x << " " << a33x << "  | " << a31inv << " " << a32inv << " " << a33inv << std::endl;
		std::cout << "[KERNEL APROX.] Matrix [" << a11x << ", " << a12x << ", " << a13x << "; " \
																																										<< a21x << ", " << a22x << ", " << a23x << "; " \
																																										<< a31x << ", " << a32x << ", " << a33x << "]" << std::endl;

		std::cout << "[KERNEL APROX.] Inv(Matrix) [" << a11inv << ", " << a12inv << ", " << a13inv << "; " \
																																															<< a21inv << ", " << a22inv << ", " << a23inv << "; " \
																																															<< a31inv << ", " << a32inv << ", " << a33inv << "]" << std::endl;

		*/



	//Compute ghost note velocity
	dens_temp = rhoWsum*a11inv + drhodx*a12inv + drhody*a13inv;

	real grx = rhoWsum*a21inv + drhodx*a22inv + drhody*a23inv;
	real gry = rhoWsum*a31inv + drhodx*a32inv + drhody*a33inv;

	dens =( dens_temp + grx*gnd.x + gry*gnd.y);

	//std::cout << " --> Case A: ACTIVE." << std::endl;

	}
	else if (a11x>0)
	{

	dens = rhoWsum/a11x;

	//// VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//// //VEL_APROX = VEL_APROX + vWsum; //??? make it better

	//std::cout << " --> Case B: ACTIVE." << std::endl;

	}
	else
	{

		dens = 1000;

	//std::cout << " --> Case C: ACTIVE." << std::endl;

	}

	return dens;


}

/* Approximate velocity with kernel aproximation to certain point */
realvec Kernel_velocity_approximation_with_mass
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	realvec VEL_APROX = {0. , 0.};

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;

	unsigned int npx, npy;
	unsigned int ncx, ncy;

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr = {0., 0.};
	realvec nv = {0., 0.};
	realvec dr = {0., 0.};
	real drs = 0;
	real nrho = 0;

	double *kernel;

	realvec vWsum = {0., 0.};
	real Wsum = 0;
	realvec vvv = {0., 0.};
	real mWsum = 0;
	//realvec vWsum;
	//real Wsum ;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nv = particles.data.v[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		vWsum = vWsum + nv*kernel[0]; // += not overloaded yet
		vvv = vvv + nv*kernel[0]*(m/nrho);
		mWsum += m*kernel[0]/nrho;

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//VEL_APROX.y = mWsum;
	//VEL_APROX = VEL_APROX + vWsum; //??? make it better
	//VEL_APROX = VEL_APROX + vvv; //??? make it better
	//VEL_APROX =  vvv; //??? make it better

	if(Wsum == 0){VEL_APROX = {0. ,0.};}

	return VEL_APROX;
}


/* Approximate velocity with kernel aproximation to certain point */
real Kernel_positionVec_divergence
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	real divr = 0.;

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	real nrho;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	real rhoWsum = 0;
	real Wsum = 0;
	real rdWsum = 0;

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += nrho*kernel[0]; // += not overloaded yet

		rdWsum += (dr.x*dr.x*kernel[1] + dr.y*dr.y*kernel[1])*m/nrho;



		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	//dens =  rhoWsum / Wsum; //??? make it better
	//if(Wsum == 0){dens = 1000;}
	divr = rdWsum;

	return divr;
}

realvec Kernel_particles_concentration_gradient
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	realvec dC = {0., 0.};

	/* Get cell and nb's cells for point of aprotixmation. */
	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	real h = particles.data_const.h;
	real m = particles.data_const.m;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	npx = (int)((ar.x - simulation_data.x_0 + of)/kh);
	npy = (int)((ar.y - simulation_data.y_0 + of)/kh);

	idx c = POS(npx, npy, ncx, ncy);
	idx zz, zp, pp, pz, pm, zm, mm, mz, mp;
	std::vector<int> ac;

	//Load indices of neighbour cells
	zz = c;
	zp = c + 1;
	pp = c + 1 + ncy;
	pz = c + ncy;
	pm = c - 1 + ncy;
	zm = c - 1;
	mm = c - 1 - ncy;
	mz = c - ncy;
	mp = c + 1 - ncy;

	ac = {zz, zp, pp, pz, pm, zm, mm, mz, mp};

		/* Debug */
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> c: " << c << std::endl;
		//std::cout << "ACCELERATION VEL. KERNEL APPROX -> zz: " << zz << " zp: " << zp << " pp: " << pp << " pz: " << pz << " pm: " << pm << " zm: " << zm << " mm: " << mm << " mz:" << mz << " mp: " << mp << std::endl;

	realvec nr;
	real nrho;
	realvec dr;
	real drs;
	//real drdv;

	double *kernel;

	real rhoWsum = 0;
	real Wsum = 0;
	//real rdWsum = 0;
	realvec dcSum = {0., 0.};

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

			/* this shloud he here I guess */
		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}

		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		nr = particles.data.r[particles.cells[cl].cp[n]];
		nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		rhoWsum += nrho*kernel[0]; // += not overloaded yet

		dcSum = dcSum + dr*kernel[1]*m/nrho;



		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	//dens =  rhoWsum / Wsum; //??? make it better
	//if(Wsum == 0){dens = 1000;}
	dC = dcSum;

	return dC;
}
