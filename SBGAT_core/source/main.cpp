
#include <iostream>
#include "ShapeModelImporter.hpp"
#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"
#include "Densities.hpp"
#include <chrono>


#include <armadillo>

int main( int argc, char** argv ) {

	std::chrono::time_point<std::chrono::system_clock> start, end;
    

	ShapeModel shape_model;
	ShapeModelImporter shape_io("../../resources/KW4Alpha.obj");
	shape_io.load_shape_model(&shape_model);

	DynamicAnalyses dynamic_analyses(&shape_model);
    start = std::chrono::system_clock::now();
	
	dynamic_analyses.compute_pgm(Density::KW4_ALPHA_DENSITY);
	

    end = std::chrono::system_clock::now();
 	std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	return 0;
}
