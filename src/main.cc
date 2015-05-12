#include <iostream>
#include "MGSolver.hh"
#include "Smoother.hh"

int 
main(int argc, char **args)
{
	if(argc != 3){
		std::cout << "Usage: ./mgsolve <number_of_levels> <number_of_V-cycles>" << std::endl;
	return 0;
	}
	
	int l;
	int n;

	std::istringstream iss(args[1]);
	if(!(iss >> l)){
		std::cerr << "Could not parse number of level argument: " << args[1] << std::endl;
		return 1;
	}
	iss.str("");
	iss.clear();

	iss.str(args[2]);
	if(!(iss >> n)){
		std::cerr << "Could not parse number of V-cycle argument: " << args[2] << std::endl;
		return 1;
	}
//	iss.str("");
//	iss.clear();

	Smoother smoother;
	MGSolver solver(l, smoother);
	solver.initialize_assignment_01();
//	solver.initialize_random();
	solver.v_cycle(2, 1, n);
	solver.saveToFile();

	return 0;
}

