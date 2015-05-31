#include <iostream>
#include "MGSolver.hh"
#include "Smoother.hh"
#include <sys/time.h>


	int 
main(int argc, char **args)
{
	if(argc != 3){
		std::cout << "Usage: ./mgsolve <number_of_levels> <number_of_V-cycles>" << std::endl;
		return 0;
	}

	int l;
	int n;

	//time
	//	siwir::Timer ti;
	//	double time;

	struct timeval t0, t;

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

	Smoother smoother;
	MGSolver solver(l, smoother);
#ifdef NEUMANN
	solver.initialize_seminar();
#elseif SEMINAR
	solver.initialize_seminar();
#else
	solver.initialize_seminar();
#endif
	solver.saveToFile("init.dat");

	gettimeofday(&t0, NULL);

	solver.v_cycle(2, 1, n);

	gettimeofday(&t, NULL);
	std::cout << "Wall clock time of MG execution: " << ((int64_t) (t.tv_sec - t0.tv_sec) * (int64_t)1000000 + 					(int64_t)t.tv_usec - (int64_t)t0.tv_usec) * 1e-3 << " ms" << std::endl;

	//	time = ti.elapsed();
	//	std::cout << "Time: " << "\t" << time << std::endl;

	solver.saveToFile("solution.dat");

	return 0;
}

