#include "SPH_defs.h"
#include "SPH_particle_system.h"

		Create_Wall_x_bt(particles, 0. - dp*2, 2. + dp*2, 0. , -1, {0., 1.});

		//DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.1 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);

		DRAW_CreateInletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {2.0, 0.}, inlet);
		DRAW_CreateOutletY(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {2.0, 0.}, outlet, true);

	//	DRAW_CreateInletY(particles, 0.9  - 4*dp, 0.0 + dp, 0.2, +1, 4, {-2.0, 0.}, inlet);
		//DRAW_CreateOutletY(particles, 0.1  + 4*dp, 0.0 + dp, 0.2, -1, 4, {-2.0, 0.}, outlet, true);

		//Create_Block_of_Fluid(particles, 0.2 + dp, 0.1 + dp, 0.8, 0.3);
