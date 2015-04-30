#ifndef MGSOLVER_H
#define MGSOLVER_H

#include <vector>

#include "Array.hh"
#include "Types.hh"
#include "Smoother.hh"

class MGSolver
{
	public:
	
		// will allocate memory for all needed arrays
		MGSolver ( int levels, Smoother & smoother );

		// initializes rhs and initial values for ex01
		void initialize_assignment_01();

		// performes times v-cycles
		void v_cycle ( int pre_smooth,
                       int post_smooth,
                       int times
                     );
	

	private:

		int levels_;
	
		std::vector<Array *> v_grids_;
		std::vector<Array *> r_grids_;
		std::vector<real>    h_intervals_;    // grid spacing for each level

		Smoother smoother_;

		void apply_operator_2d_poisson ( Array & source, Array & target );

		void v_cycle_pvt ( int pre_smooth,
                           int post_smooth,
                           int level			// solve when level == 1
                         );

		real residual_2d ( Array & u,
                           Array & f,
                           real h
                         );


};

#endif



