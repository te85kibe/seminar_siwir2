
#include "MGSolver.hh"
#include "Smoother.hh"

int 
main()
{

	Smoother smoother;
	MGSolver solver(4, smoother);
	solver.initialize_assignment_01();
//	solver.initialize_random();
	solver.v_cycle(2, 1, 10);

	return 0;
}

