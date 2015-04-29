#ifndef MGSOLVER_H
#define MGSOLVER_H

#include <vector>

#include "Array.hh"
#include "Types.hh"

class MGSolver
{
	public:
	
		// will allocate memory for all needed arrays
		MGSolver ( int levels );

		// performes times v-cycles
		void v_cycle ( int pre_smooth,
                       int post_smooth,
                       int times
                     );
	

	private:

		int levels_;
	
		std::vector<Array *> v_grids_;
		std::vector<Array *> r_grids_;

		void apply_operator_2d_poisson ( Array & source, Array & target );

		void v_cycle_pvt ( int pre_smooth,
                           int post_smooth,
                           int level			// solve when level == 1
                         );

};

#endif



