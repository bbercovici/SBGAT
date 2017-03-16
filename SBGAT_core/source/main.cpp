
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
	ShapeModelImporter shape_io("../../resources/earth_spherical.obj");
	shape_io.load_shape_model(&shape_model);

	DynamicAnalyses dynamic_analyses(&shape_model);
	
	// dynamic_analyses.compute_pgm(Density::KW4_ALPHA_DENSITY);
	arma::vec p = {6378000,0,0};
    start = std::chrono::system_clock::now();
	
	std::cout << dynamic_analyses.pgm_acceleration(p,Density::EARTH_DENSITY);

    end = std::chrono::system_clock::now();
 	std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";


	return 0;
}
