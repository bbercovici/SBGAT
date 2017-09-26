#include <iostream>
#include <armadillo>
#include <chrono>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>
#include <Constants.hpp>

int main( int argc, char** argv ) {


	SBGAT_CORE::FrameGraph frame_graph;
	SBGAT_CORE::ShapeModel shape_model("B", &frame_graph);
	SBGAT_CORE::ShapeModelImporter shape_io("../eros.obj", 1);

	shape_io.load_shape_model(&shape_model);

	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	arma::mat Cnm_total;
	arma::mat Snm_total;


	dynamic_analyses.compute_exterior_sh_coefs_normalized(
	    Cnm_total,
	    Snm_total,
	    10,
	    16,
	    2670000000000.0,
	    false);

	int n = 10;
	int m = 9;

	std::cout << "\nC" << std::to_string(n) << std::to_string(m) << " :" << Cnm_total.row(n)(m) << std::endl;;


	std::cout << "S" << std::to_string(n) << std::to_string(m) << " :" << Snm_total.row(n)(m);

	return 0;
}
