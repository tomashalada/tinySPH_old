#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_find_pairs.h"
#include "SPH_data_to_vtk.h"
#include "SPH_grid_to_vtk.h"
#include "SPH_dataInterpolation.h"
#include "SPH_pressure.h"
#include "SPH_interaction.h"
#include "SPH_integration.h"
#include "SPH_output_info.h"
#include "draw_geometry.h"
#include "SPH_mDBC.h"
#include "SPH_btdebug.h"

//#include "SPH_moving_boundary.h"
#include "SPH_inlet_outlet.h" //inlet, outlet, functions
#include "SPH_inlet_outlet_data.h" //inlet, outlet zones structures
#include "SPH_buffer_zones.h"
//#include "SPH_shifting.h"

struct SPH_simulation
{

	#include "constants.h"

	// Simulation objects and shared structures
	Particle_system particles;
	Simulation_data simulation_data;

	void INIT
	()
	{

		// Draw case
		#include "geometry.h"

		// Create grid for co-interacting pairs
		Divide_To_Cells(particles, simulation_data);
		// Initial contidion to VTK
		write_to_ASCII_VTK(particles, fileName_initCond);
		// Find neihgbours test
		Particle_To_Cell(particles, simulation_data);
		// Grid to VTK
		write_mesh_to_ASCII_VTK(particles, simulation_data, fileName_grid);
		// Simulation info
		Output_file(particles.data_const, simulation_data, fileName_info);

	}


	void RUN
	()
	{

		// Main time loop
		while(step < step_end)
		{
		step++;
		std::cout << "SIMULATION -> RUN: Step: " << step << std::endl;

		/* Apply mDBC, compute pressure */
		mDBC_compute_density_mdbcGeo(particles, simulation_data);
		// Density_to_pressure(particles);
		std::cout << "SIMULATION -> RUN: Apply mDBC, compute pressure. DONE. " << std::endl;

		// /* Integration */
		// Integrate_SymplecticPredictor(particles, dt);
		// if(step > 1 || step%40 != 0){Integrate_Verlet(particles, dt);}
		// Integrate_Euler(particles, dt);
		// std::cout << "SIMULATION -> RUN: First part of integraion. DONE. " << std::endl;

		/* Clear cells */
		Clear_cells(particles);
		std::cout << "SIMULATION -> RUN: Clear cells. DONE. " << std::endl;

		/* Particles to cells (find neihgbours) */
		Particle_To_Cell(particles, simulation_data);
		std::cout << "SIMULATION -> RUN: Particles to cells. DONE. " << std::endl;

		/* Integrate densiy and compute pressure */
		Density_to_pressure(particles);
		std::cout << "SIMULATION -> RUN: Integrate density. DONE. " << std::endl;

		/* Compute forces */
		Compute_Forces(particles);
		std::cout << "SIMULATION -> RUN: Compute forces. DONE. " << std::endl;

		/* Integration */
		if(step > 1 || step%40 != 0){Integrate_Verlet(particles, dt);}
		Integrate_Euler(particles, dt);
		std::cout << "SIMULATION -> RUN: First part of integraion. DONE. " << std::endl;

		/* Remove partiocles out of domain */
	 RemoveParticlesOutOfDomain(particles, simulation_data);
		std::cout << "SIMULATION -> RUN: Remove particles out of domain. DONE. " << std::endl;

		#pragma omp barrier

		/* Integrate inlet and outlet */
		BZ_Integrate_inlet(particles, 0, dt);
		std::cout << "SIMULATION -> RUN: Inlet integraion. DONE. " << std::endl;
		BZ_Integrate_outlet(particles, 1, dt);
		std::cout << "SIMULATION -> RUN: Outlet integraion. DONE. " << std::endl;

		/* Output to vtk. */
		if(step%save_output_interval == 0)
		{

			std::string output_file_name, output_file_namefo, output_file_name_p_out;
			output_file_name = fileName_resultsAll + std::to_string(step) + ".vtk";
			output_file_namefo = fileName_resultsFluidOnly + std::to_string(step) + ".vtk";
			output_file_name_p_out = fileName_particlesOut + std::to_string(step) + ".vtk";

			write_to_ASCII_VTK(particles, output_file_name);
			//write_to_ASCII_VTK_noIOzones(particles, output_file_namefo);
 		//excluded_particle_write_to_ASCII_VTK(particles, output_file_name_p_out);

		}


		std::cout << "SIMULATION -> STEP RECAP <-" << std::endl;
		std::cout << "np: " << particles.np << " | Data array size -> r: " << particles.data.r.size() << " rho: " << particles.data.rho.size() << std::endl;
		std::cout << "np_out: " << particles.np_out << std::endl;


		/* Interpolate and save interpolated data. */
		if(step%save_output_interval_interpolated == 0)
		{

			std::string output_file_nameInterpol;
			output_file_nameInterpol = fileName_resultsInterpol + std::to_string(step) + ".vtk";

			//void GenerateInterpol
			//(Particle_system &particles, Simulation_data simulation_data, std::string fname, int step, real x0, real y0, real xm, real ym)
			GenerateInterpol(particles, simulation_data, output_file_nameInterpol, 0.0, 0., 0.8, 0.7);
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
	: particles(hh, kap, cs, 1000., visco, delta, dp),
			simulation_data(x_0, x_1, y_0, y_1)
	{

	}

};
