
#include <iostream>

#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"
#include "Densities.hpp"



#include <armadillo>

int main( int argc, char** argv ) {

	ShapeModel shape_model;
	shape_model.load("../../resources/KW4alpha.obj");

	DynamicAnalyses dynamic_analyses(&shape_model);
	
	dynamic_analyses.compute_pgm(Density::KW4_ALPHA_DENSITY);

	return 0;
}
