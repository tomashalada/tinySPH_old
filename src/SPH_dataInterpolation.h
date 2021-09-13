#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_data.h"
#include "SPH_kernel_approx.h"

#include "SPH_grid_to_vtk.h"

void GenerateInterpol
(Particle_system &particles, Simulation_data simulation_data, std::string filename, real x0, real y0, real xm, real ym)
{

	real h = particles.data_const.h;
	//Domain limits

	/*
	 x0 = 0.1;
	 y0 = 0.0;
	 xm = 0.8;
	 ym = 0.2;
	*/

	/*
	x0 = 0.87;
	y0 = 0.0;
	xm = 1.5;
	ym = 0.2;
	*/

	real sf = 2; //Smoothness factor
	real dh = h/sf;

	int ncx = (int)((xm - x0)/h) * sf;
	int ncy = (int)((ym - y0)/h) * sf;
	std::cout << "[DATA INTERPOLATION] ncx = " << ncx << " ncy = " << ncy << std::endl;

	std::vector<realvec> vel(ncx*ncy);
	std::vector<realvec> pos(ncx*ncy);

	//realvec r = {x0, y0};


	//Create grid over data
	#pragma omp parallel for schedule(static) collapse(2)
	for(int y = 0; y < ncy; y++)
		for(int x = 0; x < ncx; x++){
	{

		realvec r;
		r.x = x0 + x*dh;
		r.y = y0 + y*dh;
		//pos[POS(x, y, ncx, ncy)] = r;
		//pos[POS(x,y,ncx,ncy)] = r;

		vel[POS(x, y, ncx, ncy)] = Kernel_velocity_approximation_with_mass(particles, simulation_data, r);
		//vel[POS(y, x, ncx, ncy)].y = Kernel_velocity_approximation(particles, simulation_data, r).x;
		//vel[POS(x, y, ncx, ncy)].x = POS(x, y, ncx, ncy);
		//vel[POS(x, y, ncx, ncy)].y = POS(y, x, ncx, ncy);

		}
	}
	#pragma omp barrier

	write_INTERPOLATED_mesh_to_ASCII_VTK(ncx, ncy, dh, x0, y0, vel, filename);
	//write_to_ASCII_VTK_POINTS("OUTPUTOO.vtk", vel, pos);
	//write_INTERPOLATED_DAT(pos, vel, "INTERPOLATED_DAT", step);

	//write_INTERPOLATED_HM(pos, vel, "INTERPOLATED_HM", step, ncx, ncy);


}
