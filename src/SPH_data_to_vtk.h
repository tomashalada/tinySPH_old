#include "SPH_particle_system.h"
#include "SPH_defs.h"

void write_to_ASCII_VTK
(Particle_system particles, std::string output_file_name)
{

	unsigned int np = particles.np;

	std::ofstream file;
	file.open(output_file_name);

	// file header
	file << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;

	// point coordinates list
	file << "POINTS " << np << " float" << std::endl;
	for(const realvec &r: particles.data.r)
	{
		file << r.x << " " << r.y << " 0" << std::endl;
	}

	//	data fields
	file << "POINT_DATA " << np << std::endl << "FIELDS FieldData 4" << std::endl;

	file << "type 1 " << np << " int" << std::endl;
	for(const idx &type: particles.data.part_type)
	{
		file << type << std::endl;
	}

	file << "pressure 1 " << np << " float" << std::endl;
	for(const real &p: particles.data.p)
	{
		file << p << std::endl;
	}

	file << "density 1 " << np << " float" << std::endl;
	for(const real &rho: particles.data.rho)
	{
		file << rho << std::endl;
	}

	file << "velocity 3 " << np << " float" << std::endl;
	for(const realvec &v : particles.data.v){
		file << v.x << " " << v.y << " 0" << std::endl;
	}

	file.close();

}

void excluded_particle_write_to_ASCII_VTK
(Particle_system particles, std::string output_file_name)
{

	unsigned int np = particles.np_out;

	std::ofstream file;
	file.open(output_file_name);

	// file header
	file << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;

	// point coordinates list
	file << "POINTS " << np << " float" << std::endl;
	for(const realvec &r: particles.particles_out.r)
	{
		file << r.x << " " << r.y << " 0" << std::endl;
	}

	//	data fields
	file << "POINT_DATA " << np << std::endl << "FIELDS FieldData 4" << std::endl;

	file << "type 1 " << np << " int" << std::endl;
	for(const idx &type: particles.particles_out.part_type)
	{
		file << type << std::endl;
	}

	file << "pressure 1 " << np << " float" << std::endl;
	for(const real &p: particles.particles_out.p)
	{
		file << p << std::endl;
	}

	file << "density 1 " << np << " float" << std::endl;
	for(const real &rho: particles.particles_out.rho)
	{
		file << rho << std::endl;
	}

	file << "velocity 3 " << np << " float" << std::endl;
	for(const realvec &v : particles.particles_out.v){
		file << v.x << " " << v.y << " 0" << std::endl;
	}

	file.close();

}


void write_to_ASCII_VTK_noIOzones
(Particle_system particles, std::string output_file_name)
{

	unsigned int np = particles.np;
	unsigned int fp = 0;

	for(int i = 0; i < np; i++)
	{
		if((particles.data.part_type[i] == fluid)|| (particles.data.part_type[i] == wall)){ fp++;	}
	}

	std::ofstream file;
	file.open(output_file_name);

	// file header
	file << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;

	// point coordinates list
	file << "POINTS " << fp << " float" << std::endl;
	for(int i = 0; i < np; i++)
	{
		if((particles.data.part_type[i] == fluid)|| (particles.data.part_type[i] == wall)){
		file << particles.data.r[i].x << " " << particles.data.r[i].y << " 0" << std::endl;
		}
	}

	//	data fields
	file << "POINT_DATA " << fp << std::endl << "FIELDS FieldData 3" << std::endl;

	file << "pressure 1 " << fp << " float" << std::endl;
	for(int i = 0; i < np; i++)
	{
		if((particles.data.part_type[i] == fluid)|| (particles.data.part_type[i] == wall)){
		file << particles.data.p[i] << std::endl;
		}
	}

	file << "density 1 " << fp << " float" << std::endl;
	for(int i = 0; i < np; i++)
	{
		if((particles.data.part_type[i] == fluid)|| (particles.data.part_type[i] == wall)){
		file << particles.data.rho[i] << std::endl;
		}
	}

	file << "velocity 3 " << fp << " float" << std::endl;
	for(int i = 0; i < np; i++)
	{
		if((particles.data.part_type[i] == fluid)|| (particles.data.part_type[i] == wall)){
		file << particles.data.v[i].x << " " << particles.data.v[i].y << " 0" << std::endl;
		}
	}

	file.close();

}

