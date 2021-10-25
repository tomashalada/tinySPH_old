#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"

void write_mesh_to_ASCII_VTK
(Particle_system particles, Simulation_data sd, std::string filename)
{

	std::ofstream file;

	int ncx, ncy;
	ncx = particles.pairs.ncx + 1;
	ncy = particles.pairs.ncy + 1;

	/* This should be in some data struct. */
	real kh, of;
	kh = particles.data_const.h*particles.data_const.kap;
	of = particles.data_const.h*0.273;
	//double of = 0;

	const int nFields = 2;
	int nCells = (ncx - 1 ) * (ncy - 1);

	file.open(filename);

	//hlavička souboru
	file << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET STRUCTURED_GRID" << std::endl
	<< "DIMENSIONS "<< ncx << " " << ncy << " 1" << std::endl << "POINTS " << ncx*ncy << " float" << std::endl;

	//vypsání bodů sítě do souboru
	for(int i=0; i<ncx; i++){
		for(int j=0; j<ncy; j++){
			file << sd.x_0 + i*kh - of << " " << sd.y_0 + j*kh - of << " 0" << std::endl;
		}
	}


	file << "CELL_DATA " << nCells << std::endl << "FIELD FieldData " << nFields << std::endl;


	//jednotlivé pole s hodnotama

	//skalární pole
	file << "NumberOfParticles 1 " << nCells << " int" << std::endl;
	for(int x=0; x<ncx-1; x++) {
		for(int y=0; y<ncy-1; y++){

		// na každém řádku hodnota buňky
			file << particles.cells[POS(x,y, ncx-1, ncy-1)].np << std::endl;


		}
	}

	//skalární pole
	file << "Cell1DIndex 1 " << nCells << " int" << std::endl;
	for(int x=0; x<ncx-1; x++) {
		for(int y=0; y<ncy-1; y++){

		// na každém řádku hodnota buňky
			file << POS(x,y, ncx-1, ncy-1) << std::endl;


		}
	}

}

void write_buffer_mesh_to_ASCII_VTK
(Particle_system particles,int ncx, int ncy, double kh, double x_0, double y_0, std::vector<int> bfc, int step)
{

	std::ofstream file;

	/*
	int ncx, ncy;
	ncx = particles.pairs.ncx + 1;
	ncy = particles.pairs.ncy + 1;
	*/

	/* This should be in some data struct. */
	//real kh, of;
	//kh = particles.data_const.h*particles.data_const.kap;
	//of = particles.data_const.h*0.273;
	//double of = 0;

	const int nFields = 2;
	int nCells = (ncx - 1 ) * (ncy - 1);

	std::string fnm = "io_mesh/outlet_mesh_";
	fnm += std::to_string(step) + ".vtk";
	file.open(fnm);
	fnm = "io_mesh/outlet_mesh_";

	//hlavička souboru
	file << "# vtk DataFile Version 2.4" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET STRUCTURED_GRID" << std::endl
	<< "DIMENSIONS "<< ncx << " " << ncy << " 1" << std::endl << "POINTS " << ncx*ncy << " float" << std::endl;

	//vypsání bodů sítě do souboru
		for(int j=0; j<ncy; j++){
	for(int i=0; i<ncx; i++){
			file << x_0 + i*kh - kh/2<< " " << y_0 + j*kh - kh/2<< " 0" << std::endl;
		}
	}


	file << "CELL_DATA " << nCells << std::endl << "FIELD FieldData " << nFields << std::endl;


	//jednotlivé pole s hodnotama

	//skalární pole
	file << "NumberOfParticles 1 " << nCells << " int" << std::endl;
		for(int y=0; y<ncy-1; y++){
	for(int x=0; x<ncx-1; x++) {

		// na každém řádku hodnota buňky
			file << bfc[POS(x,y, ncx-1, ncy-1)] << std::endl;


		}
	}

	//skalární pole
	file << "Cell1DIndex 1 " << nCells << " int" << std::endl;
	for(int x=0; x<ncx-1; x++) {
		for(int y=0; y<ncy-1; y++){

		// na každém řádku hodnota buňky
			file << POS(x,y, ncx-1, ncy-1) << std::endl;


		}
	}

}

void write_INTERPOLATED_mesh_to_ASCII_VTK
(int ncx, int ncy, double kh, double x_0, double y_0, std::vector<realvec> v, std::vector<real> rho, std::vector<real> press, std::string filename)
{

	std::ofstream file;

	/*
	int ncx, ncy;
	ncx = particles.pairs.ncx + 1;
	ncy = particles.pairs.ncy + 1;
	*/

	/* This should be in some data struct. */
	//real kh, of;
	//kh = particles.data_const.h*particles.data_const.kap;
	//of = particles.data_const.h*0.273;
	//double of = 0;
	int pocetBoduSiteX = ncx;
	int pocetBoduSiteY = ncy;

	int sitParaviewX = pocetBoduSiteX + 1;
	int sitParaviewY = pocetBoduSiteY + 1;
	int pocBunek = (pocetBoduSiteX)*(pocetBoduSiteY);
	sitParaviewX = pocetBoduSiteX;
	sitParaviewY = pocetBoduSiteY;
	const int pocPoli = 3;

	std::cout << "[DATA INTERPOLATION - GRID TO VTK] ncx = " << ncx << " ncy = " << ncy << std::endl;

	//std::string fnm = "RESULT_INTERPOLATEDhr_";
	file.open(filename);


	file << "# vtk DataFile Version 3.0" << std::endl;
	file << "vtk output" << std::endl;
	file << "ASCII" << std::endl;
	file << "DATASET STRUCTURED_GRID" << std::endl;
	file << "DIMENSIONS "<< sitParaviewX << " " << sitParaviewY << " 1" << std::endl;
	file << "POINTS " << sitParaviewX*sitParaviewY << " float" << std::endl;

		for(float j=0;j<sitParaviewY;j++)
	for(float i=0;i<sitParaviewX;i++)
			file << x_0+ i*kh << " " << y_0 + j*kh << " 0" << std::endl;

/*
Hlavička souboru .vtk, sem nesahat!
(říká, že teď budeme vypisovat hodnoty)
*/
	file << "POINT_DATA " << pocBunek << std::endl;
	file << "FIELD FieldData " << pocPoli << std::endl;



	//jednotlivé pole s hodnotama

	/*
	//skalární pole
	file << "NumberOfParticles 1 " << nCells << " int" << std::endl;
		for(int y=0; y<ncy-1; y++){
	for(int x=0; x<ncx-1; x++) {

		// na každém řádku hodnota buňky
			file << bfc[POS(x,y, ncx-1, ncy-1)] << std::endl;


		}
	}
	*/

	//skalární pole
/*
	file << "VELOCITY 3 " << nCells << " double" << std::endl;
	for(int x=0; x<ncy; x++) {
		for(int y=0; y<ncx; y++){

		// na každém řádku hodnota buňky
			file << v[POS(x,y, ncx, ncy)].x << " " << v[POS(x,y, ncx-1, ncy-1)].y << " 0" << std::endl;
			//file << x << " " << y << " 0" << std::endl;

		}
	}
	*/
	file << "DENSITY 1 " << pocBunek << " float" << std::endl;
		for(int j=0;j<pocetBoduSiteY;j++){
	for(int i=0;i<pocetBoduSiteX;i++){
			//	file << data[i*pocetBoduSiteY +j] << std::endl; //<------------------ TADY, pozor na indexování
			file << rho[POS(i,j, pocetBoduSiteX, pocetBoduSiteY)] << std::endl;
	}}

	file << "PRESSURE 1 " << pocBunek << " float" << std::endl;
		for(int j=0;j<pocetBoduSiteY;j++){
	for(int i=0;i<pocetBoduSiteX;i++){
			//	file << data[i*pocetBoduSiteY +j] << std::endl; //<------------------ TADY, pozor na indexování
			file << press[POS(i,j, pocetBoduSiteX, pocetBoduSiteY)] << std::endl;
	}}

	file << "VELOCITY 3 " << pocBunek << " float" << std::endl;
		for(int j=0;j<pocetBoduSiteY;j++){
	for(int i=0;i<pocetBoduSiteX;i++){
			//	file << data[i*pocetBoduSiteY +j] << std::endl; //<------------------ TADY, pozor na indexování
			file << v[POS(i,j, pocetBoduSiteX, pocetBoduSiteY)].x << " " << v[POS(i,j, pocetBoduSiteX, pocetBoduSiteY)].y << " 0" << std::endl;
}}

file.close();

}

//   void write_to_ASCII_VTK_POINTS
//   (std::string output_file_name, std::vector<realvec> vel, std::vector<realvec> pos)
//   {
//
//   	unsigned int np = pos.size();
//
//   	std::ofstream file;
//   	file.open(output_file_name);
//
//   	// file header
//   	file << "# vtk DataFile Version 3.0" << std::endl << "vtk output" << std::endl << "ASCII" << std::endl << "DATASET POLYDATA" << std::endl;
//
//   	// point coordinates list
//   	file << "POINTS " << np << " float" << std::endl;
//   	for(const realvec &r: pos)
//   	{
//   		file << r.x << " " << r.y << " 0" << std::endl;
//   	}
//
//   	//	data fields
//   	file << "POINT_DATA " << np << std::endl << "FIELDS FieldData 1" << std::endl;
//
//   	/*
//   	file << "type 1 " << np << " int" << std::endl;
//   	for(const idx &type: particles.data.part_type)
//   	{
//   		file << type << std::endl;
//   	}
//
//   	file << "pressure 1 " << np << " float" << std::endl;
//   	for(const real &p: particles.data.p)
//   	{
//   		file << p << std::endl;
//   	}
//
//   	file << "density 1 " << np << " float" << std::endl;
//   	for(const real &rho: particles.data.rho)
//   	{
//   		file << rho << std::endl;
//   	}
//   	*/
//
//   	file << "velocity 3 " << np << " float" << std::endl;
//   	for(const realvec &v : vel){
//   		file << v.x << " " << v.y << " 0" << std::endl;
//   	}
//
//   	file.close();
//
//   }
//
//   void write_INTERPOLATED_DAT
//   (std::vector<realvec> r, std::vector<realvec> v, std::string fname, int step)
//   {
//
//   	std::ofstream file;
//
//   	std::string fnm = fname;
//   	fnm += std::to_string(step) + ".vtk";
//   	file.open(fnm);
//
//   	for(unsigned int i = 0; i < r.size(); i++){
//   		file << r[i].x << " " << r[i].y << " " << v[i].x << " " << v[i].y << std::endl;
//   	}
//
//
//   }
//
//   void write_INTERPOLATED_HM
//   (std::vector<realvec> r, std::vector<realvec> v, std::string fname, int step, int ncx, int ncy)
//   {
//
//   	std::ofstream file;
//
//   	std::string fnm = fname;
//   	fnm += std::to_string(step) + ".vtk";
//   	file.open(fnm);
//
//   	for(unsigned int y = 0; y < ncy; y++){
//   	for(unsigned int x = 0; x < ncx; x++){
//   		file << v[POS(x,y,ncx,ncy)].x << " ";
//
//   	}
//   	file << std::endl;
//   	}
//
//
//   }
