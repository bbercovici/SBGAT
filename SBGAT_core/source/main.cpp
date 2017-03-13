
#include <iostream>

#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"
#include "Densities.hpp"
#include <chrono>


#include <armadillo>

int main( int argc, char** argv ) {

	std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    

	ShapeModel shape_model;
	shape_model.load("../../resources/KW4alpha.obj");



	DynamicAnalyses dynamic_analyses(&shape_model);
	
	dynamic_analyses.compute_pgm(Density::KW4_ALPHA_DENSITY);
    end = std::chrono::system_clock::now();
 	std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	return 0;
}
