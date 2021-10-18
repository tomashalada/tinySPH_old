#pragma once

#include "SPH_defs.h"

//double *Wendland_kernel(double r, double h){
//
//	static double smooth_fc[3];
//	double rs = fabs(r / h);
//
//	//double dim_param = 5/(14*M_PI*h*h);
//	double dim_param = 7/(4*M_PI*h*h);// * 1.7951958;
//
//	if (rs <= 2.0){
//		smooth_fc[0] = (pow((1-rs/2),4)*(2*rs + 1)) * dim_param;
//		smooth_fc[1] = -2.7852/(h*h*h)*(pow((1-rs/2),3))*rs/r;
//	}
//	else{
//		smooth_fc[0] = 0.0;
//		smooth_fc[1] = 0.0;
//	}
//	if (r == 0.0){
//		smooth_fc[0] = dim_param;
//		smooth_fc[1] = 0.0;
//	}
//	return smooth_fc;
//}


double *Wendland_kernel(double r, double h){

	static double smooth_fc[3];
	double rs = fabs(r / h);

	//double dim_param = 5/(14*M_PI*h*h);
	double dim_param = 7/(4*M_PI*h*h);// * 1.7951958;

	if (rs <= 2.0){
		smooth_fc[0] = ((1-rs/2)*(1-rs/2)*(1-rs/2)*(1-rs/2)*(2*rs + 1)) * dim_param;
		smooth_fc[1] = -2.7852/(h*h*h)*(1-rs/2)*(1-rs/2)*(1-rs/2)/r;
	}
	else{
		smooth_fc[0] = 0.0;
		smooth_fc[1] = 0.0;
	}
	if (r == 0.0){
		smooth_fc[0] = dim_param;
		smooth_fc[1] = 0.0;
	}
	return smooth_fc;
}
