#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_find_pairs.h"
#include "SPH_data_to_vtk.h"
#include "SPH_grid_to_vtk.h"

#include "SPH_density.h"
#include "SPH_pressure.h"
#include "SPH_acceleration.h"

#include "SPH_integration.h"
#include "SPH_moving_boundary.h"

#include "SPH_output_info.h"
#include "draw_geometry.h"

#include "SPH_mDBC.h"

#include "SPH_inlet_outlet.h" //inlet, outlet, functions
#include "SPH_inlet_outlet_data.h" //inlet, outlet zones structures

// template <typename TYPES>
struct SPH_simulation
{

	/*** Simulation parameters & settings ***/

	double x_0 = 0.; // temp
	double x_1 = 2.; // temp
	double y_0 = 0.; // temp
	double y_1 = 1.; // temp

	real dpLR = 0.01; //temp
	real dpHR = 0.0075; //temp
	//real dpHR = 0.005; //temp
	real dt = 0.0001; //temp
	real kap = 1.;

	unsigned int save_output_interval = 25; //temp, default 25
	unsigned int step = 0; //Time in steps
	unsigned int step_end = 5001; //Endtime in steps;

 std::string output_file_nameLR = "resultsT/particles_";
	//std::string output_file_name_p_out = "p_out/particles_";

	/*** Inlet and outlet zone settings ***/

	//Inlet constant buffer in LR region.
	realvec LRinlet_vel = {3. , -0.5};
	real LRinlet_x0 = 0.1 + 4*dpLR;
	real LRinlet_y0 = 0.0 + dpLR;

	//Outlet dynamic buffer in LR region, with dynamic height
	real LRoutletDyn_x0 = 0.9-4*dpLR;
	real LRoutletDyn_y0 = 0.0 + dpLR;
	real LRoutletDyn_xm = 0.9;
	real LRoutletDyn_ym = 0.2;

	/*** Help & temp variables ***/

	double hhLR = 1*sqrt(2*dpLR*dpLR);
	double hhHR = 1*sqrt(2*dpLR*dpHR);
	//double hh = 0.02;
	int nvlLR = std::ceil(kap*hhLR/dpLR);
	int nvlHR = std::ceil(kap*hhHR/dpHR);

	/*** Simulation objects and shared structures ***/
	Particle_system particlesLR;
	Particle_system particlesHR;

	Simulation_data simulation_dataLR;
	Simulation_data simulation_dataHR;

	//outletcheat
	double cheat = 0;
	double htz = 0;

	void INIT
	()
	{

		/* Draw case */
		//Draw_geometry_dam_break(particles, dp);
		//Draw_geometry_piston(particles, dp);
		//Draw_geometry_dam_break_with_inlet(particles, dp);
		//Draw_geometry_piston_narrow(particles, dp);
		Draw_geometry_dam_break_with_inlet_and_outletLRSR(particlesLR, dpLR);

		// Create grid for co-interacting pairs
		Divide_To_Cells(particlesLR, simulation_dataLR);

		// Initial contidion to VTK
		//write_to_ASCII_VTK(particlesLR, "init_conditionLR.vtk");

		/* Test */
		Particle_To_Cell(particlesLR, simulation_dataLR);

		// Grid to VTK
		//write_mesh_to_ASCII_VTK(particlesHR, simulation_dataLR);

		// Simulation info
		//Output_file(particlesLR.data_const, simulation_dataLR);

	}


	void RUN
	()
	{

		/* Main time loop. */
		while(step < step_end)
		{

		step++;
		std::cout << "SIMULATION -> RUN: Step: " << step << std::endl;

		/* First part of integration */
		if(step > 1)
		{

			Integrate_LeapFrog_partOne(particlesLR, dt);

		}
		std::cout << "SIMULATION -> RUN: First part integraion. DONE. " << std::endl;

		/* Moving boundary, TEMP */
		// Moving_boundary(particlesLR, dpLR); // remove dp
		// Moving_boundary(particlesHR, dpHR); // remove dp
		// std::cout << "SIMULATION -> RUN: Moving boundary. DONE. " << std::endl;

		/* Clear cells, TEMP */
		Clear_cells(particlesLR);
		std::cout << "SIMULATION -> RUN: Clear cells. DONE. " << std::endl;

		/* Particles to cells (find neihgbours) */
		Particle_To_Cell(particlesLR, simulation_dataLR);
		std::cout << "SIMULATION -> RUN: Particles to cells. DONE. " << std::endl;

		/* Compute density change */
		Compute_Density(particlesLR);
		mDBC_compute_density(particlesLR, simulation_dataLR, dpLR/2);
		std::cout << "SIMULATION -> RUN: Compute density. DONE. " << std::endl;

		/* Integrate densiy and compute pressure, CHANGE */
		Integrate_density_compute_pressure(particlesLR, dt);
		std::cout << "SIMULATION -> RUN: Integrate density. DONE. " << std::endl;

		/* Compute acceleration */
		Compute_Acceleration(particlesLR);
		std::cout << "SIMULATION -> RUN: Compute acceleration. DONE. " << std::endl;

		/* Second part of integration */
		Integrate_LeapFrog_partTwo(particlesLR, dt, step);
		std::cout << "SIMULATION -> RUN: Second part integraion. DONE. " << std::endl;

		/*** Update inlet buffer ***/
		//void Integrate_inlet(Particle_system &particles, double dt, double x_b, realvec v_inl, real b_len)
		//Integrate_inlet(particlesLR, dt, LRinlet_x0, LRinlet_vel, 4*dpLR);
		Integrate_inlet_wide(particlesLR, dt, LRinlet_x0, 0.2, LRinlet_vel, 4*dpLR, 0.2);

		//std::cout << "Debug: " << particlesLR.zones.size() << " Zone coordinates: [x0, y0] = [" << particlesLR.zones[0].x0 << "," << particlesLR.zones[0].y0 << "], [xm, ym] = [" << particlesLR.zones[0].xm << "," << particlesLR.zones[0].ym << "]"<< ", bl: " << particlesLR.zones[0].bl <<std::endl;

		//Integrate_inlet(particlesLR, particlesLR.zones[0], dt);
		std::cout << "SIMULATION -> RUN: Update inlet LR buffer. DONE. " << std::endl;

		#pragma omp barrier
								/* Update outlet buffer */
								//Outlet_treatment_update_velocity(particles, simulation_data, 0.9-4*dp -dp/2, dp); // -dp/2 nDBC mirror line

								/* Static buffer  - WRONG*/
								//Outlet_treatment_remove_fluid(particles, 0.9-4*dp, 4*dp, dt, dp); // this works somehow

								/* Dynamic buffer - RIGHT */

//		//void Outlet_treatment_remove_fluid_with_renew(Particle_system &particles, double x_b, double b_len, double dt, double dp)
//		Outlet_treatment_remove_fluid_with_renew(particlesLR, LRoutletDyn_x0, 4*dpLR, dt, dpLR); // this works somehow
//
//		//vodi Renew_buffer(Particle_system &particles, Simulation_data simulation_data, double x0_b, double y0_b, double xm_b, double ym_b, double dp, int step, double &cheat)
//		Renew_buffer(particlesLR, simulation_dataLR, LRoutletDyn_x0, LRinlet_y0, LRoutletDyn_xm, LRoutletDyn_ym, dpLR, step, cheat);
//
//		Outlet_treatment_update_velocity(particlesLR, simulation_dataLR, 0.9-4*dpLR -dpLR/2, dpLR); // -dp/2 nDBC mirror line
//									//Outlet_treatment_update_velocity(particlesLR, simulation_dataHR, 0.9-4*dpHR -dpHR/2, dpHR); // -dp/2 nDBC mirror line
//
//
//																			/* Interpolate inlet velocityies*/
//																			//Renew_transition_zone(particlesHR, cheat, htz, 0.9, 0.0 + dpHR, -1);
//																			Renew_transition_zone_smooth(particlesHR, cheat, htz, 0.9, 0.0 + dpHR, -1);
//
//																			//fromDraw//////////////////////Inlet_y_direction(particles, 0.9, 0.0 + dp, 0.2, -1, 4, {0.0, 0.});
//																			Treat_transition_inlet(particlesHR, particlesLR, simulation_dataHR, simulation_dataLR);
//
//																			/////////Integrate_inlet(particlesHR, dt, 0.9 - 4*dpHR, {2.0, 0.}, 4*dpHR); // buffer possitions???? Fixit
//																			Integrate_inlet(particlesHR, dt, 0.9, {0.0, 0.}, 4*dpHR); // buffer possitions???? Fixit
//																			std::cout << "SIMULATION -> RUN: Second part integraion. DONE. " << std::endl;


		// --- Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});
		std::cout << "SIMULATION -> RUN: Inlet and outlet treatment. DONE. " << std::endl;

		/* Output to vtk. */
		if(step%save_output_interval == 0)
		{

			output_file_nameLR += std::to_string(step) + ".vtk";

			//output_file_name_p_out += std::to_string(step) + ".vtk";

			write_to_ASCII_VTK(particlesLR, output_file_nameLR);

 		//excluded_particle_write_to_ASCII_VTK(particles, output_file_name_p_out);

		}

		//Clear output file name variable for next output step
 	output_file_nameLR = "resultsT/particles_"; //move
		//output_file_name_p_out = "p_out/particles_";


		std::cout << "SIMULATION -> STEP RECAP <-" << std::endl;
		std::cout << "np: " << particlesLR.np << " | Data array size -> r: " << particlesLR.data.r.size() << " rho: " << particlesLR.data.rho.size() << std::endl;
		std::cout << "np_out: " << particlesLR.np_out << std::endl;

		} // Main time loop.



	} // RUN function

	// Compute and insert sim. data to proper structures
	void PREP_SIMULATION_DATA
	()
	{

		//LR
		simulation_dataLR.hh = hhLR;
		simulation_dataLR.kh = hhLR*kap;
		simulation_dataLR.nvl = nvlLR;

		simulation_dataLR.x_0 -= nvlLR*dpLR;
		simulation_dataLR.x_m += nvlLR*dpLR;

		simulation_dataLR.y_0 -= nvlLR*dpLR;
		simulation_dataLR.y_m += nvlLR*dpLR;

		//HR
		simulation_dataHR.hh = hhHR;
		simulation_dataHR.kh = hhHR*kap;
		simulation_dataHR.nvl = nvlHR;

		simulation_dataHR.x_0 -= nvlHR*dpHR;
		simulation_dataHR.x_m += nvlHR*dpHR;

		simulation_dataHR.y_0 -= nvlHR*dpHR;
		simulation_dataHR.y_m += nvlHR*dpHR;

	}

	/* Init */
	SPH_simulation
	()
	: particlesLR(hhLR, kap, 44.29, 1000., dpLR),
	 	particlesHR(hhHR, kap, 44.29, 1000., dpHR),
			simulation_dataLR(x_0, x_1, y_0, y_1),
			simulation_dataHR(x_0, x_1, y_0, y_1)
	{

	/* Debug */
	std::cout << "SIMULATION -> Smoothing length h_dp: " << hhHR << std::endl;


	}

};
