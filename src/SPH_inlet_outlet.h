#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_kernel_approx.h"

#include "SPH_grid_to_vtk.h"
#include "SPH_inlet_outlet_data.h"


/* Inlet treatment */
void Integrate_inlet
(Particle_system &particles, double dt, double x_b, realvec v_inl, real b_len)
{

		//#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] != inlet){continue;}
			//particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt;
			particles.data.r[p].x = particles.data.r[p].x + particles.data.v[p].x * dt;

			if(particles.data.r[p].x > x_b)
			{

				particles.data.part_type[p] = fluid;
				particles.Add_particle(inlet, 0, 1000, {particles.data.r[p].x - b_len, particles.data.r[p].y}, v_inl, {0., 0.});

			}

		}

} // function


void Integrate_inlet_wide
(Particle_system &particles, double dt, double x_b, double y_b, realvec v_inl, real b_len, real b_height)
{

		//#pragma omp parallel for schedule(static)
		for(int p = 0; p < particles.np; p++)
		{

			if(particles.data.part_type[p] != inlet){continue;}
			particles.data.r[p].x = particles.data.r[p].x + particles.data.v[p].x * dt;

			if(particles.data.r[p].x > x_b)
			{

				particles.data.part_type[p] = fluid;
				particles.Add_particle(inlet, 0, 1000, {particles.data.r[p].x - b_len, particles.data.r[p].y}, v_inl, {0., 0.});

			}

		}

} // function





void Outlet_treatment_update_velocity
(Particle_system &particles, Simulation_data simulation_data, double x_b, double dp)
{

	realvec r;
	realvec gn; //ghost node coordinates
	realvec gnd; //ghost node distance
	realvec vbp; //velocity of buffer particle
	real rhobp = 0;

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

/* Make grid over OUTLET (wokrs with inlet also) buffer zone
			and  check if every cell with length of edge dp contains particle.
			If not so, create one buffer particle in center of empty cell. */
void Renew_buffer
(Particle_system &particles, Simulation_data simulation_data, double x0_b, double y0_b, double xm_b, double ym_b, double dp, int step, double &cheat)
{

	real bh = 0;
	real hh = particles.data_const.h;
	real kap = particles.data_const.kap;

	//real bhc = 0; real bhl = 0; real bhr = 0;
	//bhc = Water_Elevation(particles, simulation_data, 0.9-4*dp - 1.5*hh , 0 + hh*kap, dp);
	//bhl = Water_Elevation(particles, simulation_data, 0.9-4*dp - 2*hh , 0 + hh*kap, dp);
	//bhr = Water_Elevation(particles, simulation_data, 0.9-4*dp - 1*hh , 0 + hh*kap, dp);
	//bh = std::max({bhc, bhl, bhr});

	bh = Water_Elevation(particles, simulation_data, 0.9-4*dp - 2*hh , 0 + hh*kap, dp);

	/* Debug */
	std::cout << "[RENEW BUFFER] Near buffer watter elevation: " << bh  << std::endl;
	//std::ofstream fileWL;
	//fileWL.open("waterLevel.dat", std::ios_base::app);
	//fileWL << step << " " << bh << " " << bhl << " " << bhc << " " << bhr << std::endl;
	//fileWL.close();

	if(bh > cheat){cheat = bh;} //TEMP
	else{bh = cheat;} //TEMP

	real lenx = xm_b - x0_b;
	//real leny = ym_b - y0_b; //fixed outlet height
	real leny = bh - y0_b; //variable buffer length

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


/* Retype fluid to outlet if particle crosses buffer region,
			update outlet particles (integration) and remove outlet
			particles if crosses oposit direction of buffer. */
void Outlet_treatment_remove_fluid_with_renew
(Particle_system &particles, double x_b, double b_len, double dt, double dp)
{

	//#pragma omp parallel for schedule(static)
	for(int p = 0; p < particles.np; p++)
	{

		if(particles.data.part_type[p] == fluid)
		{

			if(particles.data.r[p].x > x_b) // change fluid particle to buffer particle
			{

				particles.data.part_type[p] = outlet; //retype fluid incomming to buffer
				particles.data.a[p] = {0., 0};
				particles.data.drho[p] = 0.;

			}

		} // particle fluid if
		else if(particles.data.part_type[p] == outlet)
		{

			particles.data.r[p] = particles.data.r[p] + particles.data.v[p] * dt; //update outlet p. position

			if(particles.data.r[p].x > x_b + b_len || particles.data.r[p].y < dp/2) // change fluid particle to buffer particle
			{

				particles.Remove_particle(p); //remove particle outgoing from buffer

			}

		} /// particle outlet if

	} // particle loop

} // function

void Treat_transition_inlet
(Particle_system &particlesHR, Particle_system &particlesLR, Simulation_data simulation_dataHR, Simulation_data simulation_dataLR)
{

	real rhoip;
	realvec vip;

	for(int p = 0; p < particlesHR.np; p++)
	{

		if(particlesHR.data.part_type[p] != inlet){continue;}

		vip = Kernel_velocity_approximationTRZ(particlesLR, simulation_dataLR, particlesHR.data.r[p]);
		rhoip = Kernel_density_approximationTRZ(particlesLR, simulation_dataLR, particlesHR.data.r[p]);
		////vip.x = Kernel_velocity_approximationTRZ(particlesLR, simulation_dataLR, particlesHR.data.r[p]).x;
		//std::cout << "TRANSITION ZONE -> vip.x: " << vip.x << std::endl;
		//std::cout <<  "TRANSITION ZONE inlet particle r: [" << particlesHR.data.r[p].x << "," << particlesHR.data.r[p].y << "]" << std::endl;
		////vip.y = 0;
		particlesHR.data.v[p] = vip;
		particlesHR.data.rho[p] = rhoip;

	}

}

void Renew_transition_zone
(Particle_system &particles, double &cheat, double &htz, double b_x0 ,double b_y0, int orientation)
{

	real dp = particles.data_const.dp;
	int nil_cheat = 0;
	int nil_htz = 0;
	realvec r;
	//real htz = htzz - 0.005;

	//if(cheat > 0 && htz == 0){

	//}
	//else if (cheat > 0 && htz >= 0){
	std::cout << "OUTPUT ZONE: cheat: " << cheat << " htz: " << htz << std::endl;
	if (cheat > 0 && htz >= 0){
			nil_cheat = (int)(cheat/dp);
			nil_htz = (int)(htz/dp);
			//nil_cheat = std::round(cheat/dp);
			//nil_htz = std::round(htz/dp);

			std::cout << "OUTPUT ZONE: Active - cheat: " << cheat << " htz: " << htz << std::endl;
			std::cout << "OUTPUT ZONE: Active - nil_cheat: " << nil_cheat << " nil_htz: " << nil_htz << std::endl;
			if (nil_cheat > nil_htz){

																					for(int i = 0; i < (nil_cheat - nil_htz); i++)
																					{
																					r.y = b_y0 + nil_htz*dp + (i )*dp ;

																					// 	// 	//real dp = particles.data_const.dp;
																					// 	// 	//real kh = particles.data_const.h * particles.data_const.kap;
																					// 	// 	//realvec r, v, a;

																					// 	// 	////gets length of wall line and number of particles for the wall
																					// 	// 	//real len = LEN(x, y_0, x, y_1);
																					// 	// 	//int npart = (int)(len/dp);
																					// 	// 	///*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

																					// 	// 	//// number of virtual layers
																					// 	// 	//int nvl = std::ceil(kh/dp);
																					// 	// 	////nvl = 1;


																					//for(int l = 0; l < nvl; l++)
																					for(int l = 0; l < 4; l++)
																					{

																						r.x = b_x0 + dp*l*orientation;


																							//particles.Add_particle(inlet, 0, particles.data_const.rho0, r, v, a);
																							particles.Add_particle(inlet, 0, 1000, r, {0., 0.}, {0., 0.});
																					}
																				}

			}
			htz = cheat;

	}

}


void Renew_transition_zone_smooth
(Particle_system &particles, double &cheat, double &htz, double b_x0 ,double b_y0, int orientation)
{

	real dp = particles.data_const.dp;
	int nil_cheat = 0;
	int nil_htz = 0;
	realvec r;
	double bh_dif; //buffer height difference
	int ntg; //number of layers to generate

	//if(cheat > 0 && htz == 0){

	//}
	//else if (cheat > 0 && htz >= 0){
	std::cout << "OUTPUT ZONE: cheat: " << cheat << " htz: " << htz << std::endl;
	if (cheat > 0 && htz >= 0){
			nil_cheat = (int)(cheat/dp) ;
			nil_htz = (int)(htz/dp) ;
			//nil_cheat = std::round(cheat/dp);
			//nil_htz = std::round(htz/dp);

			std::cout << "OUTPUT ZONE: Active - cheat: " << cheat << " htz: " << htz << std::endl;
			std::cout << "OUTPUT ZONE: Active - nil_cheat: " << nil_cheat << " nil_htz: " << nil_htz << std::endl;
			if (nil_cheat > nil_htz){

																					//bh_dif = cheat - htz;
																					//ntg = std::ceil(bh_dif/dp);

																					for(int i = 0; i < (nil_cheat - nil_htz); i++)
																					//for(int i = 0; i < ntg; i++)
																					{
																					r.y = b_y0 + nil_htz*dp + (i )*dp ;

																					// 	// 	//real dp = particles.data_const.dp;
																					// 	// 	//real kh = particles.data_const.h * particles.data_const.kap;
																					// 	// 	//realvec r, v, a;

																					// 	// 	////gets length of wall line and number of particles for the wall
																					// 	// 	//real len = LEN(x, y_0, x, y_1);
																					// 	// 	//int npart = (int)(len/dp);
																					// 	// 	///*if(modulo != 0){Delka hrany nekoresponduje se zadanym NP}*/

																					// 	// 	//// number of virtual layers
																					// 	// 	//int nvl = std::ceil(kh/dp);
																					// 	// 	////nvl = 1;


																					//for(int l = 0; l < nvl; l++)
																					for(int l = 0; l < 4; l++)
																					{

																						r.x = b_x0 + dp*l*orientation;


																							//particles.Add_particle(inlet, 0, particles.data_const.rho0, r, v, a);
																							particles.Add_particle(inlet, 0, 1000, r, {0., 0.}, {0., 0.});
																					}
																				}

			}
			htz = cheat;

	}

}

