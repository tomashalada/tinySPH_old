#pragma once

#include "SPH_defs.h"

enum bt
{

	inlet_v_const = 0,
	inlet_v_profile = 1,
	inlet_v_height = 2,

	outlet_v_const = 3,
	outlet_v_profile = 4,
	outlet_v_height = 5

};

struct BufferZone
{

	realvec orientation; //{+1,0},{-1,0},{0,+1},{0,-1}
	real x0, y0, xm, ym; //zone region
	unsigned int nbl; //number of buffer layers

	real bl, bh; //buffer length, buffer
	unsigned int btype; //buffer type
	realvec v_b; //default buffer velocity

	/* Help, test, special */
	real bh_old; //buffer height in previous time step

	BufferZone
	(real _x0, real _xm, real _y0, real _ym, idx nbl, realvec _v0, realvec _orientation, unsigned int _btype) :
	x0(_x0), xm(_xm), y0(_y0), ym(_ym), v_b(_v0), orientation(_orientation), btype(_btype)
	{

		bl = fabs(xm - x0);
		bh = fabs(ym - y0);

		if(orientation.x != 0){ bl = fabs(xm - x0); bh = fabs(ym - y0); }
		else if(orientation.y != 0){ bl = fabs(ym - y0); bh = fabs(xm - x0); }
		else{ std::cout << "[CREATE BUFFER ZONE] >> INVALID INLET ORIENTATION." << std::endl; } //throw error

	}

};




//struct with information about buffer zones
struct BufferZones
{

	real x0, y0, xm, ym;//zone region
	real bl = xm - x0;//buffer lenght
	real bh = ym - y0;//buffer height
	unsigned int nbl;//number of buffer layers

	realvec v_b; //initial buffer velocity

	BufferZones
	(real _x0, real _xm, real _y0, real _ym, idx nbl, realvec _v0) :
	x0(_x0), xm(_xm), y0(_y0), ym(_ym), v_b(_v0) {}

};

struct InletBuffer : BufferZones
{

	unsigned short int part_type = inlet;
	using BufferZones::BufferZones;

};

struct OutletBuffer : BufferZones
{

	unsigned short int part_type = outlet;
	using BufferZones::BufferZones;

};
