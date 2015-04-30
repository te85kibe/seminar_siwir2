
#include "MGSolver.hh"
#include "Smoother.hh"

int 
main()
{

	Smoother smoother;
	MGSolver solver(4, smoother);
	solver.initialize_assignment_01();
	solver.v_cycle(1, 1, 100);

	return 0;
}

