#pragma once

#include "SPH_defs.h"
#include "SPH_particle_system.h"
#include "SPH_create.h"

//#include "gen_case.h" //dimensions and basic parameters

void Draw_geometry_dam_break
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 1. + dp*2, 0., -1);
		Create_Wall_x(particles, 0. - dp*2, 1. + dp*2, 1., +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

		Create_Wall_y(particles, 0., 0 + dp, 1.0 - dp, -1);
		Create_Wall_y(particles, 1., 0 + dp, 1.0 - dp, +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);

}

void Draw_geometry_piston
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		/* Vnejsi stena */
		Create_Wall_x(particles, 0. - dp*2, 1.0 + dp*2, 0., -1);
		Create_Wall_y(particles, 0., 0 + dp, 1.0 - dp, -1);
		Create_Wall_y(particles, 1.0, 0 + dp, 1.0 - dp, +1);

		/* Vnitrni stena */
		Create_Wall_x(particles, 0.2 + dp*2, 0.8 - dp*2, 0.2, +1); // vnitrni stena == opacna orientace
		Create_Wall_y(particles, 0.2, 0.2, 1.0, +1); // vnitrni stena == opacna orientace
		Create_Wall_y(particles, 0.8, 0.2, 1.0, -1); // vnitrni stena == opacna orientace

		Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.2 - dp, 0.7);
		Create_Block_of_Fluid(particles, 0.2, 0 + dp, 0.8, 0.2 - dp);
		Create_Block_of_Fluid(particles, 0.8 + dp, 0 + dp, 1.0, 0.3);

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);

}

void Draw_geometry_piston_narrow
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		/* Vnejsi stena */
		Create_Wall_x(particles, 0. - dp*2, 0.5 + dp*2, 0., -1);
		Create_Wall_y(particles, 0., 0 + dp, 1.0 - dp, -1);
		Create_Wall_y(particles, 0.5, 0 + dp, 1.0 - dp, +1);

		/* Vnitrni stena */
		Create_Wall_x(particles, 0.2 + dp*2, 0.3 - dp*1, 0.2, +1); // vnitrni stena == opacna orientace
		Create_Wall_y(particles, 0.2, 0.2, 1.0, +1); // vnitrni stena == opacna orientace
		Create_Wall_y(particles, 0.3, 0.2, 1.0, -1); // vnitrni stena == opacna orientace

		Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.2 - dp, 0.7);
		Create_Block_of_Fluid(particles, 0.2, 0 + dp, 0.3 + dp, 0.2 - dp);
		Create_Block_of_Fluid(particles, 0.3 + dp, 0 + dp, 0.5, 0.3);

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);

}

void Draw_geometry_dam_break_with_inlet
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 1. + dp*2, 0., -1);
		Create_Wall_x(particles, 0. - dp*2, 1. + dp*2, 1., +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

		Create_Wall_y(particles, 0., 0 + dp, 1.0 - dp, -1);
		Create_Wall_y(particles, 1., 0 + dp, 1.0 - dp, +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);
		Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2 + dp, +1, 4, {2.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outlet
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

		Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
		Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);
		Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		//Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outletLR
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
						//disabled for mDBC Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

						//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
						//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);
		//Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outletHR
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
							//disabled Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

							//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
							//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);
		//////////////Inlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		///////////////////////Inlet_y_direction(particles, 0.9, 0.0 + dp, 0.2, -1, 4, {0.0, 0.});
		//Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outletLR_PF
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
						//disabled for mDBC Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		Create_Wall_x(particles, 0., 2., 0.2, +1);

						//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
						//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2 -dp, -1, 4, {2.0, 0.}, inlet);
		Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2 -dp, +1, 4, {0.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outletHR_PF
(Particle_system &particles, real dp, real dpLR)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
							//disabled Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		Create_Wall_x(particles, 0., 2., 0.2, +1);

							//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
							//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Create moving wall
		//Create_Wall_x_mov(particles, 0.3 + dp, 0 + dp, 0.6, +1);
		Inlet_y_direction(particles, 0.9  - 4*dpLR, 0.0 + dp, 0.2 - dp, -1, 4, {0.0, 0.});
		//------DRAW_CreateInletY(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, -1, 4, {0.0, 0.}, inlet);
		///////////////////////Inlet_y_direction(particles, 0.9, 0.0 + dp, 0.2, -1, 4, {0.0, 0.});
		//Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});

}

void Draw_geometry_dam_break_with_inlet_and_outletLRSR
(Particle_system &particles, real dp)
{

		// SEQUENCE TO CREATE PARTICLES
		Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0. , -1);
						//disabled for mDBC Create_Wall_x(particles, 0. - dp*2, 2. + dp*2, 0.4 + dp, +1);
		//Create_Wall_x(particles, 0., 1., 1., +1);

						//disabled Create_Wall_y(particles, 0., 0 + dp, 0.4 , -1);
						//disabled Create_Wall_y(particles, 2., 0 + dp, 0.4 , +1);

		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 1-dp, 0.5);
		//Create_Block_of_Fluid(particles, 0 + dp, 0 + dp, 0.3, 0.5);
		//Create_Block_of_Fluid_with_speed(particles, 0.15, 0 + dp , 0.85, 0.2, {2. , 0.});

		//Inlet_y_direction(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.});
		//ORIGO DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);
		DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.2 + dp, 0.4, -1, 4, {3.0, -0.5}, inlet);
		//Outlet_y_direction(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.});

}
