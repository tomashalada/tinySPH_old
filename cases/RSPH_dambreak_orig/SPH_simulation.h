#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_find_pairs.h"
#include "SPH_data_to_vtk.h"
#include "SPH_grid_to_vtk.h"
#include "SPH_dataInterpolation.h"

#include "RSPH_density.h"
#include "SPH_density.h"
#include "SPH_pressure.h"
#include "RSPH_acceleration.h"
#include "SPH_acceleration.h"

//#include "SPH_integration.h"
#include "RSPH_integration.h"
#include "SPH_moving_boundary.h"

#include "SPH_output_info.h"
//#include "draw_geometry.h"
#include "SPH_create_bt.h"

#include "SPH_mDBC.h"

#include "SPH_inlet_outlet.h" //inlet, outlet, functions
#include "SPH_inlet_outlet_data.h" //inlet, outlet zones structures

#include "SPH_shifting.h"
#include "SPH_btdebug.h"

struct SPH_simulation
{

	#include "constants.h"

	/*** Simulation objects and shared structures ***/
	Particle_system particles;

	Simulation_data simulation_data;

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
		//Draw_geometry_dam_break_with_inlet_and_outletLR(particles, dp);
		#include "geometry.h"

		// Create grid for co-interacting pairs
		Divide_To_Cells(particles, simulation_data);

		// Initial contidion to VTK
		write_to_ASCII_VTK(particles, fileName_initCond);

		/* Test */
		Particle_To_Cell(particles, simulation_data);

		// Grid to VTK
		write_mesh_to_ASCII_VTK(particles, simulation_data);

		// Simulation info
		Output_file(particles.data_const, simulation_data, fileName_info);

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

			Integrate_LeapFrog_partOne(particles, dt);

		}
		std::cout << "SIMULATION -> RUN: First part integraion. DONE. " << std::endl;

		/* Moving boundary, TEMP */
		// Moving_boundary(particles, dp);
		// Moving_boundary(particlesHR, dpHR);
		// std::cout << "SIMULATION -> RUN: Moving boundary. DONE. " << std::endl;

		/* Clear cells, TEMP */
		Clear_cells(particles);
		std::cout << "SIMULATION -> RUN: Clear cells. DONE. " << std::endl;

		/* Particles to cells (find neihgbours) */
		Particle_To_Cell(particles, simulation_data);
		std::cout << "SIMULATION -> RUN: Particles to cells. DONE. " << std::endl;

		/* Compute density change */
		RSPHCompute_Density(particles);
		//Compute_Density(particles);
		//mDBC_compute_density(particles, simulation_data, dp/2);
		std::cout << "SIMULATION -> RUN: Compute density. DONE. " << std::endl;

		/* Integrate densiy and compute pressure, CHANGE */
		Integrate_density_compute_pressure(particles, dt);
		std::cout << "SIMULATION -> RUN: Integrate density. DONE. " << std::endl;

		/* Compute acceleration */
		RSPHCompute_Acceleration(particles);
		//Compute_Acceleration(particles);
		std::cout << "SIMULATION -> RUN: Compute acceleration. DONE. " << std::endl;

		/* Second part of integration */
		Integrate_LeapFrog_partTwo(particles, dt, step);
		std::cout << "SIMULATION -> RUN: Second part integraion. DONE. " << std::endl;

	 RemoveParticlesOutOfDomain(particles, simulation_data);
		std::cout << "SIMULATION -> RUN: Remove particles out of domain. DONE. " << std::endl;

		//real csm;
		//real vnorm;
		//csm = NORM(particles.data.v[0].x, particles.data.v[0].y);
		//for(int i = 1; i<particles.np; i++)
		//{

		//vnorm = NORM(particles.data.v[i].x, particles.data.v[i].y);
		//if(vnorm > csm){csm = vnorm;}

		//}
		//particles.data_const.cs = csm*10;


		// --- Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});
		std::cout << "SIMULATION -> RUN: Inlet and outlet treatment. DONE. " << std::endl;

		/* Output to vtk. */
		if(step%save_output_interval == 0)
		{

			std::string output_file_name;

			output_file_name = fileName_resultsAll + std::to_string(step) + ".vtk";

			output_file_namefo = fileName_resultsFluidOnly + std::to_string(step) + ".vtk";

			//output_file_name_p_out += std::to_string(step) + ".vtk";

			//write_to_ASCII_VTK(particles, output_file_name);
			//write_to_ASCII_VTK(particlesHR, output_file_nameHR);

			write_to_ASCII_VTK_noIOzones(particles, output_file_namefo);

 		//excluded_particle_write_to_ASCII_VTK(particles, output_file_name_p_out);

		}



		std::cout << "SIMULATION -> STEP RECAP <-" << std::endl;
		std::cout << "np: " << particles.np << " | Data array size -> r: " << particles.data.r.size() << " rho: " << particles.data.rho.size() << std::endl;
		std::cout << "np_out: " << particles.np_out << std::endl;


		//OUTPUT INTERPOLATE DATA
		if(step%save_output_interval_interpolated == 0)
		{

			std::string output_file_nameInterpol;

			output_file_nameInterpol = fileName_resultsInterpol + std::to_string(step) + ".vtk";

			//void GenerateInterpol
			//(Particle_system &particles, Simulation_data simulation_data, std::string fname, int step, real x0, real y0, real xm, real ym)
			GenerateInterpol(particles, simulation_data, output_file_nameInterpol, 0., 0., 0.81, 0.85);
			std::cout << "[INTERPOLATION - DONE and SAVED.]" << std::endl;
		}



		} // Main time loop.



	} // RUN function

	// Compute and insert sim. data to proper structures
	void PREP_SIMULATION_DATA
	()
	{

		//
		simulation_data.hh = hh;
		simulation_data.kh = hh*kap;
		simulation_data.nvl = nvl;

		simulation_data.x_0 -= nvl*dp;
		simulation_data.x_m += nvl*dp;

		simulation_data.y_0 -= nvl*dp;
		simulation_data.y_m += nvl*dp;


	}

	/* Init */
	SPH_simulation
	()
	//: particles(hh, kap, 44.29, 1000., dp),
	: particles(hh, kap, cs, 1000., visco, dp),
			simulation_data(x_0, x_1, y_0, y_1)
	{

	}

};
