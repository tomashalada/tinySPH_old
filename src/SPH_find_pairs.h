#pragma once

#include "SPH_particle_system.h"

void Divide_To_Cells
//(Particle_system particles)
(Particle_system &particles, Simulation_data simulation_data)
{

	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	//kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;
	//of = 0;

	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = (int)((simulation_data.x_m - simulation_data.x_0 + of)/kh) + 1;
	ncy = (int)((simulation_data.y_m - simulation_data.y_0 + of)/kh) + 1;

	particles.pairs.ncx = ncx;
	particles.pairs.ncy = ncy;


	for(int x = 0; x < ncx; x++)
	{
		for(int y = 0; y < ncy; y++)
		{

			//emplace_back takes the argmunets and create object in vector, faster then push_back with construcotr
			particles.cells.emplace_back(POS(x, y, ncx, ncy));

		}
	}

	/* Debug */
	std::cout << "FIND PAIRS -> Created ncx: " << ncx << " and ncy: " << ncy << std::endl;
	std::cout << "FIND PAIRS -> Domain x_0: " << simulation_data.x_0 << " x_m: " << simulation_data.x_m << std::endl;
	std::cout << "FIND PAIRS -> Domain y_0: " << simulation_data.y_0 << " y_m: " << simulation_data.y_m << std::endl;

}

/* Here some HPC magic is needed. */
void Particle_To_Cell
//(Particle_system particles)
(Particle_system &particles, Simulation_data simulation_data)
{

	real kh, of;
	kh = particles.data_const.kap*particles.data_const.h;
	of = particles.data_const.h*0.273;

	int npx, npy;
	int ncx, ncy;
	unsigned int ic; //cell index

	ncx = particles.pairs.ncx;
	ncy = particles.pairs.ncy;

	for(int p = 0; p < particles.np; p++)
	{

		npx = (int)((particles.data.r[p].x - simulation_data.x_0 + of)/kh);
		npy = (int)((particles.data.r[p].y - simulation_data.y_0 + of)/kh);

		//particles.cells[POS(npx, npy, ncx, ncy)].cp.push_back(particles.data.index[p]);
		particles.cells[POS(npx, npy, ncx, ncy)].cp.push_back(p);
		particles.cells[POS(npx, npy, ncx, ncy)].np++;

	}
}

void Clear_cells
(Particle_system &particles)
{

		#pragma omp parallel for schedule(static) collapse(2)
		for(int x = 0; x < particles.pairs.ncx; x++)
		{

			for(int y = 0; y < particles.pairs.ncy; y++)
			{

				particles.cells[POS(x, y, particles.pairs.ncx, particles.pairs.ncy)].np = 0;
				particles.cells[POS(x, y, particles.pairs.ncx, particles.pairs.ncy)].cp.clear();

			} // clear, for y loop

		} // clear, for x loop

}
