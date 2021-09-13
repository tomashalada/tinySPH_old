#pragma once

#include <stdio.h>
#include <cstdlib>
#include <math.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <omp.h>
#include <vector>
#include <limits>
#include <chrono>
#include <stdarg.h>
#include <sys/stat.h>
#include <sys/wait.h>

#include <algorithm>
#include <initializer_list>

// ###
#ifdef __CUDACC__
	#define CUDA_HOSTDEV __host__ __device__
	#define CUDA_HOSTDEV_NOINLINE CUDA_HOSTDEV __noinline__
#else
	#define CUDA_HOSTDEV
	#define CUDA_HOSTDEV_NOINLINE
#endif

// ###
#define SQ(x) ((x) * (x)) // square function; replaces SQ(x) by ((x) * (x)) in the code
#define NORM(x, y) sqrt(SQ(x) + SQ(y))
#define LEN(x0, y0, x1, y1) sqrt(SQ(x1-x0) + SQ(y1 - y0))

#define POS(x,y,X,Y) (y+(x)*(Y))
//#define POS(x,y,X,Y) (x+(y)*(X))

#define MAX( a , b) (((a)>(b))?(a):(b))
#define MIN( a , b) (((a)<(b))?(a):(b))


/* Type of particles */
enum{
	fluid=0,
	wall=1,
	virt=2,
	mov_a=3,
	mov_b=4,
	mov_c=5,
	inlet=6,
	outlet=7,
	outletf=8
};

// ###
#define PI  3.1415926535897932384
#define eps 0.01

// ###
typedef struct{
  double x,y;
}tdouble2;

inline tdouble2 TDouble2(double v){ tdouble2 p={v,v}; return(p); }
inline tdouble2 TDouble2(double x,double y){ tdouble2 p={x,y}; return(p); }
inline bool operator ==(const tdouble2& a, const tdouble2& b){ return(a.x==b.x&&a.y==b.y); }
inline bool operator !=(const tdouble2& a, const tdouble2& b){ return(a.x!=b.x||a.y!=b.y); }
inline tdouble2 operator +(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x+b.x,a.y+b.y)); }
inline tdouble2 operator -(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x-b.x,a.y-b.y)); }
inline tdouble2 operator *(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x*b.x,a.y*b.y)); }
inline tdouble2 operator /(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x/b.x,a.y/b.y)); }
inline tdouble2 operator +(const tdouble2& a, const double& b){ return(TDouble2(a.x+b,a.y+b)); }
inline tdouble2 operator -(const tdouble2& a, const double& b){ return(TDouble2(a.x-b,a.y-b)); }
inline tdouble2 operator *(const tdouble2& a, const double& b){ return(TDouble2(a.x*b,a.y*b)); }
inline tdouble2 operator /(const tdouble2& a, const double& b){ return(TDouble2(a.x/b,a.y/b)); }
inline tdouble2 MinValues (const tdouble2& a, const tdouble2& b){ return(TDouble2((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y))); }
inline tdouble2 MaxValues (const tdouble2& a, const tdouble2& b){ return(TDouble2((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y))); }

using realvec = tdouble2;
using real = double;
using idx = int;
