#pragma once

#include "SPH_defs.h"
#include "SPH_data.h"
//#include "SPH_simulation.h"

//Basic info to output file
void Output_file
(SPH_constants constants, Simulation_data simulation_data, std::string info_file_name)
{

	std::ofstream info_file;
	//std::string info_file_name = "results/simulation_info.out"; //not even close
	info_file.open(info_file_name);

	info_file << "--- SIMULATION INFO --- " << std::endl;
	info_file << "\n";
	info_file << "-> SPH constants: " << std::endl;
	info_file << "Alpha: " << constants.avisc << std::endl;
	info_file << "Cs: " << constants.cs << std::endl;
	info_file << "Dp: " << constants.dp << std::endl;
	info_file << "External forces: [" << constants.graviy.x << "," << constants.graviy.y << "]" << std::endl;
	info_file << "h: " << constants.h << std::endl;
	info_file << "Kap: " << constants.kap << std::endl;
	info_file << "m: " << constants.m << std::endl;
	info_file << "rho0: " << constants.rho0 << std::endl;

	info_file << "\n";
	info_file << "-> Simulation info: " << std::endl;
	info_file << "h = h(dp): " << simulation_data.hh << std::endl;
	info_file << "kh: " << simulation_data.kh << std::endl;
	info_file << "Number of virual layers (for walls) nvl: " << simulation_data.nvl << std::endl;

	info_file << "\n";
	info_file << "-> Domain size (with respect to cells): " << std::endl;
	info_file << "[x_min, y_min]: [" << simulation_data.x_0 << " , " <<  simulation_data.y_0 << "] "<< std::endl;
	info_file << "[x_max, y_max]: [" << simulation_data.x_m << " , " <<  simulation_data.y_m << "] "<< std::endl;

	info_file.close();

	std::cout << "Info file generated." << std::endl;

}


//Add specific info to output file
void Add_to_output_file
(std::string input)
{

}

//Write main info to output file
void Write_basic_simulation_info
()
{

}

void Write_water_level
(real step, real wl, std::string output_file_name_wl)
{
	std::ofstream fileWL;
	fileWL.open(output_file_name_wl, std::ios_base::app);
	fileWL << step << " " << wl << std::endl;
	fileWL.close();
}

void Write_pressure
(real step, real p, std::string output_file_name_p)
{
	std::ofstream fileWL;
	fileWL.open(output_file_name_p, std::ios_base::app);
	fileWL << step << " " << p << std::endl;
	fileWL.close();
}


