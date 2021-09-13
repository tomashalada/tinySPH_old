#pragma once

#include "SPH_defs.h"

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
