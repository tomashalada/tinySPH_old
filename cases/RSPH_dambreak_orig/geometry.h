#include "SPH_defs.h"
#include "SPH_particle_system.h"

		Create_Wall_x_bt(particles, 0.+dp , 0.8 - dp , 0. , -1, {0., 1.});
		Create_Wall_x_bt(particles, 0.+dp , 0.8 - dp , 0. - dp , -1, {0., 1.});
		Create_Wall_x_bt(particles, 0.+dp , 0.8 - dp , 0. - 2*dp , -1, {0., 1.});

		Create_Wall_y_bt(particles, 0. , 0. + dp, 1. , -1, {1., 0.});
		Create_Wall_y_bt(particles, 0.-dp , 0. + dp, 1. , -1, {1., 0.});
		Create_Wall_y_bt(particles, 0.-2*dp , 0. + dp, 1. , -1, {1., 0.});


		Create_Wall_y_bt(particles, 0.8 , 0. + dp , 1. , +1, {-1. , 0.});
		Create_Wall_y_bt(particles, 0.8+dp , 0. + dp , 1. , +1, {-1. , 0.});
		Create_Wall_y_bt(particles, 0.8+2*dp , 0. + dp , 1. , +1, {-1. , 0.});


		/*
		Create_corner_bt(particles, -1, -1, {0. -dp, 0.}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0.}, {sqrt(2), sqrt(2)});

		Create_corner_bt(particles, -1, -1, {0. , 0. -dp}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. , 0. -dp*2}, {sqrt(2), sqrt(2)});

		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0. -dp}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp, 0. -dp*2}, {sqrt(2), sqrt(2)});

		Create_corner_bt(particles, -1, -1, {0., 0.}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp, 0. -dp}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0. -dp*2}, {sqrt(2), sqrt(2)});
		*/

		Create_corner_bt(particles, -1, -1, {0. -dp, 0.}, {2/sqrt(5), 1/sqrt(5)});
		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0.}, {3/sqrt(10), 1/sqrt(10)});

		Create_corner_bt(particles, -1, -1, {0. , 0. -dp}, {1/sqrt(5), 2/sqrt(5)});
		Create_corner_bt(particles, -1, -1, {0. , 0. -dp*2}, {1/sqrt(10), 3/sqrt(10)});

		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0. -dp}, {3/sqrt(13), 2/sqrt(13)});
		Create_corner_bt(particles, -1, -1, {0. -dp, 0. -dp*2}, {2/sqrt(13), 3/sqrt(13)});

		Create_corner_bt(particles, -1, -1, {0., 0.}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp, 0. -dp}, {sqrt(2), sqrt(2)});
		Create_corner_bt(particles, -1, -1, {0. -dp*2, 0. -dp*2}, {sqrt(2), sqrt(2)});

		/*SECOND CORNER */

		Create_corner_bt(particles, -1, -1, {0.8 + dp, 0.}, {-2/sqrt(5), 1/sqrt(5)});
		Create_corner_bt(particles, -1, -1, {0.8 + dp*2, 0.}, {-3/sqrt(10), 1/sqrt(10)});

		Create_corner_bt(particles, -1, -1, {0.8 , 0. -dp}, {-1/sqrt(5), 2/sqrt(5)});
		Create_corner_bt(particles, -1, -1, {0.8 , 0. -dp*2}, {-1/sqrt(10), 3/sqrt(10)});

		Create_corner_bt(particles, -1, -1, {0.8 +dp*2, 0. -dp}, {-3/sqrt(13), 2/sqrt(13)});
		Create_corner_bt(particles, -1, -1, {0.8 +dp, 0. -dp*2}, {-2/sqrt(13), 3/sqrt(13)});

		Create_corner_bt(particles, +1, -1, {0.8, 0.}, {- sqrt(2), sqrt(2)});
		Create_corner_bt(particles, +1, -1, {0.8 + dp, 0. - dp}, {- sqrt(2), sqrt(2)});
		Create_corner_bt(particles, +1, -1, {0.8 + dp*2, 0. -dp*2}, {- sqrt(2), sqrt(2)});


		Create_Block_of_Fluid(particles, 0.0 + dp, 0 + dp, 0.2, 0.5);

