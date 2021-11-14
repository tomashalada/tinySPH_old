#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

/* Mass sumation to get water level */
real Sum_kernel_function
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

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


	double *kernel;

	real Wsum = 0;
	real mWsum = 0;

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
		realvec nr = particles.data.r[particles.cells[cl].cp[n]];
		real nrho = particles.data.rho[particles.cells[cl].cp[n]];

		//Position and velocity difference
		realvec dr = ar - nr;
		real drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		mWsum += m*kernel[0]/nrho;

		} // cycle over particles in neighbour cells
	} // cycle over neighbour cells

	//if(Wsum == 0){mWsum = 0;}
	return mWsum;
}

//GET WATER LEVEL
//x - position of measurement
//y0 - starting depth
real Water_Elevation
(Particle_system &particles, Simulation_data simulation_data, real x, real y0)
{
	real WL;

	const real dp = particles.data_const.dp;
	const real dp_wl = dp*0.00005; //sampling step

	realvec ar = {x, y0};
	//std::vector<real> mWsum_y;
	real mWsum_prev;

	for(int i = 0; i < int(simulation_data.y_m/dp_wl); i++)
	{
		ar.y += i*dp_wl;
		real mWsum = Sum_kernel_function(particles, simulation_data, ar);

		//if(mWsum < 0.5 && mWsum_prev < 0.5)
		//std::cout << ar.y/0.3 << " ";
		if(mWsum <= 0.5)
		{

			WL = ar.y - i*dp_wl;
			break;
		}
	}
	std::cout << " " << std::endl;


	return WL;
}


real Water_Elevation_inverse
(Particle_system &particles, Simulation_data simulation_data, real x, real y0)
{
	real WL;

	const real dp = particles.data_const.dp;
	const real dp_wl = dp*0.00005; //sampling step

	realvec ar = {x, y0};
	//std::vector<real> mWsum_y;
	real mWsum_prev;

	for(int i = 0; i < int(simulation_data.y_m/dp_wl); i++)
	{
		ar.y -= i*dp_wl;
		real mWsum = Sum_kernel_function(particles, simulation_data, ar);

		//if(mWsum < 0.5 && mWsum_prev < 0.5)
		//std::cout << ar.y/0.3 << " ";
		if(mWsum >= 0.5)
		{

			WL = ar.y + i*dp_wl;
			break;
		}
	}
	std::cout << " " << std::endl;


	return WL;
}
