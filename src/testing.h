#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel.h"

std::array<real, 4> Kernel_velocity_approximation_TEST
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	realvec VEL_APROX = {0. , 0.};
	real RHO_APROX = 0;
	real P_APROX = 0;
	std::array<real, 4> InterpolatedFields;

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
	real np = 0;

	double *kernel;

	realvec vWsum = {0., 0.};
	real rhoWsum = 0;
	real pWsum = 0;

	real Wsum = 0;
	realvec vvv = {0., 0.};
	real mWsum = 0;
	real mcWsum = 0;
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
		np = particles.data.p[particles.cells[cl].cp[n]];

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt(pow(dr.x, 2) + pow(dr.y, 2));
		//drdv = dr.x*dv.x + dr.y*dv.y;

		/* get kernel values
		double *Wendland_kernel(double r, double h) */
		kernel = Wendland_kernel(drs, h);

		Wsum += kernel[0];
		vWsum = vWsum + nv*kernel[0]; // += not overloaded yet
		rhoWsum = rhoWsum + nrho*kernel[0]; // += not overloaded yet
		pWsum = pWsum + np*kernel[0]; // += not overloaded yet
		vvv = vvv + nv*kernel[0]*(m/nrho);
		mWsum += m*kernel[0]/nrho;
		mcWsum += m*kernel[0];

		/* Debug */
		//std::cout << "KERNEL. APROX: nv: [" << nv.x << "," << nv.y << "]" << std::endl;


		} // cycle over particles in neighbour cells

	} // cycle over neighbour cells

	//VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	//VEL_APROX.x = mWsum;
	//VEL_APROX.y = mcWsum;
	if(mWsum > 0.5)
	{

	VEL_APROX = VEL_APROX + vWsum / Wsum; //??? make it better
	RHO_APROX = RHO_APROX + rhoWsum / Wsum;
	P_APROX = P_APROX + pWsum / Wsum;

	}
	else
	{
		VEL_APROX = {0., 0.};
		RHO_APROX = 0;
		P_APROX = 0;
	}

	//VEL_APROX.y = mWsum;
	//VEL_APROX = VEL_APROX + vWsum; //??? make it better
	//VEL_APROX = VEL_APROX + vvv; //??? make it better
	//VEL_APROX =  vvv; //??? make it better

	if(Wsum == 0){VEL_APROX = {0. ,0.};}

	//return VEL_APROX;
	InterpolatedFields[0] = RHO_APROX;
	InterpolatedFields[1] = P_APROX;
	InterpolatedFields[2] = VEL_APROX.x;
	InterpolatedFields[3] = VEL_APROX.y;

	return InterpolatedFields;
}


std::array<real, 4> Kernel_approximate_flow_variable
(Particle_system &particles, Simulation_data simulation_data, realvec ar)
{

	//Sum of variables to point ar
	realvec VEL_APPROX = {0. , 0.};
	real RHO_APPROX = 0;
	real P_APPROX = 0;
	std::array<real, 4> Variables_APPROX;

	//Temp variables
	real W_sum = 0;
	real rhoW_sum = 0;
	real pW_sum = 0;
	realvec vW_sum = {0., 0.};

	// Get cell and nb's cells for point of aprotixmation.
	const real kh = particles.data_const.kap*particles.data_const.h;
	const real of = particles.data_const.h*0.273;
	const real h = particles.data_const.h;
	const real m = particles.data_const.m;

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

	for(int &cl: ac)
	{

		if( cl < 0 ){continue;}
		//experiment for(int n = 0; n < particles.cells[cl].np; n++)
		for(int n = 0; n < particles.cells[cl].cp.size(); n++)
		{

		//if(particles.data.part_type[particles.cells[cl].cp[n]] != fluid){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == wall){continue;}
		//if(particles.data.part_type[particles.cells[cl].cp[n]] == inlet){continue;}

		//Load data of neighbour particle
		idx ni = particles.cells[cl].cp[n]; //actual particle index

		//Load data of neighbour particle
		const realvec nr = particles.data.r[ni]; //position of ative particle
		const realvec nv = particles.data.v[ni];; //position of ative particle
		const real nrho = particles.data.rho[ni];; //actual particle densiy
		const real np = particles.data.p[ni];

		//Interation variables
		realvec dr; //position
		real drs; //dr size

		//Position and velocity difference
		dr = ar - nr;
		drs = sqrt((dr.x*dr.x) + (dr.y*dr.y));

		//Assign kernel values
		kernel = Wendland_kernel(drs, h);

		//Add to sum
		W_sum += kernel[0];
		rhoW_sum += nrho*kernel[0];
		pW_sum += + np*kernel[0];
		vW_sum.x += nv.x*kernel[0];
		vW_sum.y += nv.y*kernel[0];

		} // cycle over particles in neighbour cells
	} // cycle over neighbour cells

	//In case of empty neighbour return zero values, otherwise CSPH
	if(W_sum == 0)
	{
		RHO_APPROX = 0.;
		P_APPROX = 0.;
		VEL_APPROX = {0. ,0.};
	}
	else
	{
		RHO_APPROX += RHO_APPROX + rhoW_sum / W_sum;
		P_APPROX +=  pW_sum / W_sum;
		VEL_APPROX.x +=  vW_sum.x / W_sum;
		VEL_APPROX.y +=  vW_sum.y / W_sum;
	}

	//Return approximated values
	Variables_APPROX[0] = RHO_APPROX;
	Variables_APPROX[1] = P_APPROX;
	Variables_APPROX[2] = VEL_APPROX.x;
	Variables_APPROX[3] = VEL_APPROX.y;

	return Variables_APPROX;
}

