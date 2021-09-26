#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel_approx.h"

#include "SPH_grid_to_vtk.h"
#include "SPH_inlet_outlet_data.h"


/* Inlet treatment */
void BZ_Integrate_inlet
(Particle_system &particles, unsigned int zoneIndex, double dt)
{

		BufferZone BUFFER = particles.zones[zoneIndex];

		real x_b;
		int orientation;
		if(BUFFER.orientation.x != 0){ x_b = BUFFER.x0; orientation = BUFFER.orientation.x; }
		else if(BUFFER.orientation.y != 0){ x_b = BUFFER.y0; orientation = BUFFER.orientation.y; }
		else{ std::cout << "[BZ_Integrate_inlet] >> INVALID INLET ORIENTATION." << std::endl; } //throw error

		real b_len = BUFFER.bl;
		realvec v_inl = BUFFER.v_b;
		unsigned int inletNumber = BUFFER.btype;

		std::cout << "[BZ_Integrate_inlet] >> Buffer zone data loaded. x_b: " << x_b << " b_len: " << b_len << " v_inl: [" << v_inl.x << "," << v_inl.y << "]" << std::endl;

		if(BUFFER.orientation.x != 0)
		{
			if(orientation < 0)
			{

				//#pragma omp parallel for schedule(static)
				for(int p = 0; p < particles.np; p++)
				{
					if(particles.data.part_type[p] != inletNumber){continue;}
					//particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;
					particles.data.r[p].x = particles.data.r[p].x + particles.data.v[p].x * dt;
					//if(particles.data.r[p].x > x_b)
					if(particles.data.r[p].x > x_b)
					{
						particles.data.part_type[p] = fluid;
						particles.Add_particle(inletNumber, 0, 1000, {particles.data.r[p].x - b_len, particles.data.r[p].y}, v_inl, {0., 0.});
					}
				}

			}//orientation if
			else if(orientation > 0)
			{

				//#pragma omp parallel for schedule(static)
				for(int p = 0; p < particles.np; p++)
				{
					if(particles.data.part_type[p] != inletNumber){continue;}
					//particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;
					particles.data.r[p].x = particles.data.r[p].x + particles.data.v[p].x * dt;
					if(particles.data.r[p].x < x_b)
					{
						particles.data.part_type[p] = fluid;
						particles.Add_particle(inletNumber, 0, 1000, {particles.data.r[p].x + b_len, particles.data.r[p].y}, v_inl, {0., 0.});
					}
				}

			}//orientation if
		}
		else if(BUFFER.orientation.y != 0)
		{
			if(orientation < 0)
			{

				//#pragma omp parallel for schedule(static)
				for(int p = 0; p < particles.np; p++)
				{
					if(particles.data.part_type[p] != inletNumber){continue;}
					//particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;
					particles.data.r[p].y = particles.data.r[p].y + particles.data.v[p].y * dt;
					//if(particles.data.r[p].x > x_b)
					if(particles.data.r[p].y > x_b)
					{
						particles.data.part_type[p] = fluid;
						particles.Add_particle(inletNumber, 0, 1000, {particles.data.r[p].x, particles.data.r[p].y - b_len}, v_inl, {0., 0.});
					}
				}

			}//orientation if
			else if(orientation > 0)
			{

				//#pragma omp parallel for schedule(static)
				for(int p = 0; p < particles.np; p++)
				{
					if(particles.data.part_type[p] != inletNumber){continue;}
					//particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;
					particles.data.r[p].y = particles.data.r[p].y + particles.data.v[p].y * dt;
					if(particles.data.r[p].y < x_b)
					{
						particles.data.part_type[p] = fluid;
						particles.Add_particle(inletNumber, 0, 1000, {particles.data.r[p].x, particles.data.r[p].y + b_len}, v_inl, {0., 0.});
					}
				}

			}//orientation if
		}
		else
		{
			std::cout << "INVALID INLET ORIENTATION." << std::endl;
		} //throw error

} // function


/* Outlet treatment */
void BZ_Integrate_outlet
//(Particle_system &particles, double x_b, double b_len, double dt, double dp)
(Particle_system &particles, unsigned int zoneIndex, double dt)
{

		//load constants
		real dp = particles.data_const.dp;

		//load buffer data
		BufferZone BUFFER = particles.zones[zoneIndex];

		real x_b, y0;
		int orientation;
		if(BUFFER.orientation.x != 0){ x_b = BUFFER.x0; orientation = BUFFER.orientation.x; }
		else if(BUFFER.orientation.y != 0){ x_b = BUFFER.y0; orientation = BUFFER.orientation.y; }
		else{ std::cout << "[BZ_Integrate_inlet] >> INVALID INLET ORIENTATION." << std::endl; } //throw error

		real b_len = BUFFER.bl;
		realvec v_inl = BUFFER.v_b;
		unsigned int outletNumber = BUFFER.btype;

		std::cout << "[BZ_Integrate_inlet] >> Buffer zone data loaded. x_b: " << x_b << " b_len: " << b_len << " v_inl: [" << v_inl.x << "," << v_inl.y << "]" << std::endl;


	//#pragma omp parallel for schedule(static)
	if(BUFFER.orientation.x != 0)
	{
		if(orientation > 0)
		{

			for(int p = 0; p < particles.np; p++)
			{
				if(particles.data.part_type[p] == fluid)
				{
					if(particles.data.r[p].x > x_b) // change fluid particle to buffer particle
					{
						particles.data.part_type[p] = outletNumber; //retype fluid incomming to buffer
						particles.data.a[p] = {0., 0};
						particles.data.drho[p] = 0.;
					}
				} // particle fluid if
				else if(particles.data.part_type[p] == outletNumber)
				{
					particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt; //update outlet p. position
					if(particles.data.r[p].x > x_b + b_len || particles.data.r[p].y < dp/2) // change fluid particle to buffer particle
					{
						particles.Remove_particle(p); //remove particle outgoing from buffer
					}
				} /// particle outlet if
			} // particle loop

		}
		else if(orientation < 0)
		{

			for(int p = 0; p < particles.np; p++)
			{
				if(particles.data.part_type[p] == fluid)
				{
					if(particles.data.r[p].x < x_b) // change fluid particle to buffer particle
					{
						particles.data.part_type[p] = outletNumber; //retype fluid incomming to buffer
						particles.data.a[p] = {0., 0};
						particles.data.drho[p] = 0.;
					}
				} // particle fluid if
				else if(particles.data.part_type[p] == outletNumber)
				{
					particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt; //update outlet p. position
					if(particles.data.r[p].x < x_b + b_len || particles.data.r[p].y < dp/2) // change fluid particle to buffer particle
					{
						particles.Remove_particle(p); //remove particle outgoing from buffer
					}
				} /// particle outlet if
			} // particle loop

		}
	}
	else if(BUFFER.orientation.y != 0)
	{
		if(orientation < 0)
		{



		}
		else if(orientation > 0)
		{



		}
	}

} // function

void BZ_Renew_outlet_buffer
//(Particle_system &particles, Simulation_data simulation_data, double x0_b, double y0_b, double xm_b, double ym_b, double dp, int step, double &cheat)
(Particle_system &particles, Simulation_data simulation_data, unsigned int zoneIndex, double dt)
{

	real abh = 0;
	real hh = particles.data_const.h;
	real kap = particles.data_const.kap;
	real dp = particles.data_const.dp;

	//real bhc = 0; real bhl = 0; real bhr = 0;
	//bhc = Water_Elevation(particles, simulation_data, 0.9-4*dp - 1.5*hh , 0 + hh*kap, dp);
	//bhl = Water_Elevation(particles, simulation_data, 0.9-4*dp - 2*hh , 0 + hh*kap, dp);
	//bhr = Water_Elevation(particles, simulation_data, 0.9-4*dp - 1*hh , 0 + hh*kap, dp);
	//bh = std::max({bhc, bhl, bhr});

		BufferZone OUTLET = particles.zones[zoneIndex];

		real x0_b, xm_b, y0_b, ym_b;
		if(OUTLET.orientation.x != 0){ x0_b = OUTLET.x0; }
		else if(OUTLET.orientation.y != 0){ y0_b = OUTLET.y0; }
		else{ std::cout << "INVALID OUTLET ORIENTATION." << std::endl; } //throw error

		real b_len = OUTLET.bl;
		real pbh = OUTLET.bh;
		realvec v_inl = OUTLET.v_b;

		abh = Water_Elevation(particles, simulation_data, 0.9-4*dp - 2*hh , 0 + hh*kap, dp);

	/* Debug */
	std::cout << "[RENEW BUFFER] Near buffer watter elevation: " << abh  << std::endl;
	//std::ofstream fileWL;
	//fileWL.open("waterLevel.dat", std::ios_base::app);
	//fileWL << step << " " << bh << " " << bhl << " " << bhc << " " << bhr << std::endl;
	//fileWL.close();

	if(abh > pbh){particles.zones[zoneIndex].bh = abh;} //TEMP
	else{abh = pbh;} //TEMP

	real lenx = xm_b - x0_b;
	//real leny = ym_b - y0_b; //fixed outlet height
	real leny = abh - y0_b; //variable buffer length

	if(lenx <= 0 || leny <= 0){return;}
	//if(lenx <= 0 || leny <= dp){return;}

	/* Debug */
	std::cout << "[RENEW BUFFER] Active... " << std::endl;

	idx nbcx = (int)((lenx)/(dp*1.5)) + 1;
	idx nbcy = (int)((leny)/(dp*1.0)) ;
	std::vector<idx> bfc(nbcx*nbcy, 0); // buffercells

	real plenx, pleny;
	idx pnx, pny;

	/* Debug */
	std::cout << "[RENEW BUFFER] Cells total: " << nbcx*nbcy << " nbcx: " << nbcx <<  " nbcy: " << nbcy << " lenx: " << lenx << " leny: " << leny << std::endl;

	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] != outlet){continue;}
		plenx = particles.data.r[p].x - (x0_b - dp/2);
		pleny = particles.data.r[p].y - (y0_b - dp/2);
		pnx = (int)(plenx/(dp*1.5));
		pny = (int)(pleny/(dp*1.0));
		bfc[POS(pnx, pny, nbcx, nbcy)]++;
		//std::cout << "Paticle added to idx: " << POS(pnx, pny, nbcx, nbcy) << " idx nump: " <<  bfc[POS(pnx, pny, nbcx, nbcy)] << " pnx: "<< pnx << " pny: "<<  pny << " plnex: "<< plenx <<" plney: "<< pleny << " r.x: " << particles.data.r[p].x << " r.y: "<< particles.data.r[p].y << std:: endl;


	} // particle loop

	//write_buffer_mesh_to_ASCII_VTK(Particle_system particles,int ncx, int ncy, double kh, double x_0, double y_0)
	//if(step % 25 == 0)
	//write_buffer_mesh_to_ASCII_VTK(particles, nbcx+1, nbcy+1, dp*1.0, x0_b, y0_b, bfc, step);
	realvec vbp = {0. , 0.};
	realvec gn = {0. , 0.};
	realvec gnd = {0. , 0.};
	realvec r = {0. , 0.};

	for(int x = 0 ; x < nbcx - 1; x++)
	{

		for(int y = 0 ; y < nbcy; y++)
		{

			std::cout << " [" << bfc[POS(x, y, nbcx, nbcy)] << "] ";
			if(bfc[POS(x, y, nbcx, nbcy)] == 0)
			{

				std::cout << "RENEW: " << std::endl;

				//particles.Add_particle(outlet, 0, 1000, {x0_b + x*dp, y0_b + y*dp}, {2.0, 0.}, {0., 0.});


			}
			else
			{
				//std::cout << "RENEWOFF: " << std::endl;

			}

		}

	}

		std::fill(bfc.begin(), bfc.end(), 0);

}

void BZ_Extrapolate_variables
//(Particle_system &particles, Simulation_data simulation_data, double x_b, double dp)
(Particle_system &particles, Simulation_data simulation_data, unsigned int zoneIndex)
{

	realvec r;
	realvec gn; //ghost node coordinates
	realvec gnd; //ghost node distance
	realvec vbp; //velocity of buffer particle
	real rhobp = 0;

	BufferZone OUTLET = particles.zones[zoneIndex];

	real x_b;
	if(OUTLET.orientation.x != 0){ x_b = OUTLET.x0; }
	else if(OUTLET.orientation.y != 0){ x_b = OUTLET.y0; }
	else{ std::cout << "INVALID OUTLET ORIENTATION." << std::endl; } //throw error

	real b_len = OUTLET.bl;
	realvec v_inl = OUTLET.v_b;

	std::cout << "[BZ_Extrapolate_variables] >> Buffer zone data loaded. x_b: " << x_b << " b_len: " << b_len << " v_inl: [" << v_inl.x << "," << v_inl.y << "]" << std::endl;

	real dp = particles.data_const.dp;

	//#pragma omp parallel for schedule(static) ///I THINK THIS CAN BE ON!
	for(int p = 0; p < particles.np; p++)
	{


		if(particles.data.part_type[p] != outlet){continue;}

		r = particles.data.r[p];
		gn.x = x_b - (particles.data.r[p].x - x_b);
		gn.y = particles.data.r[p].y;
		//gnd = particles.data.r[p] - gn;
		gnd = gn - particles.data.r[p];

		//gnd.y -= eps;
		/* Debug */
		std::cout << "OUTLET TREATMEN -> x_b: " << x_b << " gn.x: " << gn.x << " gn.y: " << gn.y << " r.x: " << particles.data.r[p].x << " r.y: " << particles.data.r[p].y << std::endl;


						// /* Approximate velocity with kernel aproximation to certain point */
						// realvec Kernel_velocity_approximation_MATRIX
						// (Particle_system &particles, Simulation_data simulation_data, realvec ar, realvec gnd)
						std::cout << "[ OUTLET TREATMENT ] Update velocity: Compute vbp for particle r = [" << r.x << "," << r.y << "]" << std::endl;
						vbp = Kernel_velocity_approximation_MATRIX(particles, simulation_data, gn, gnd);
						std::cout << " DONE. Vpx2 = [" << vbp.x << "," << vbp.y << "], gnd = [" << gnd.x << "," << gnd.y << "]" << std::endl;


		//vbp = Kernel_velocity_approximation(particles, simulation_data, gn);
		std::cout << " DONE. Vpx0 = [" << vbp.x << "," << vbp.y << "], gnd = [" << gnd.x << "," << gnd.y << "]" << std::endl;
		////vbp.x = Kernel_velocity_approximation(particles, simulation_data, gn).x;
		////vbp.y = 0;
		rhobp = Kernel_density_approximation_MATRIX(particles, simulation_data, gn, gnd);
		std::cout << " DONE. Dens matrix aprox. rhobp2 = " << rhobp << std::endl;
		rhobp = Kernel_density_approximation(particles, simulation_data, gn);
		std::cout << " DONE. Dens zeroor aprox. rhobp0 = " << rhobp << std::endl;
		//vbp = {2.0, 0.};
		particles.data.v[p] = vbp;
		//particles.data.rho[p] = 1000.;
		particles.data.rho[p] = rhobp;


	} // particle loop

} // function


