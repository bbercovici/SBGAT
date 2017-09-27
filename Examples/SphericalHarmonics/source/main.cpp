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
	SBGAT_CORE::ShapeModelImporter shape_io("../eros_64.obj", 1);

	shape_io.load_shape_model(&shape_model);


	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	arma::mat Cnm_total;
	arma::mat Snm_total;

	int degree = 5;

	dynamic_analyses.compute_exterior_sh_coefs_normalized(
		Cnm_total,
		Snm_total,
		degree,
		16,
		2670000000000.0,
		true);


	std::cout << std::setprecision(10);

	arma::mat coefs = arma::zeros<arma::mat>((degree + 1) * (degree + 2)/2 - 1,4);
	
	std::cout << Cnm_total << std::endl;


	for (int n = 1 ; n < degree + 1 ; ++n){
		for (int m = 0 ; m <= n; ++m){

			int index = (n) * (n + 1 )/2 + m - 1;

			coefs.row(index)(0) = n;
			coefs.row(index)(1) = m;
			coefs.row(index)(2) = Cnm_total.row(n)(m);
			coefs.row(index)(3) = Snm_total.row(n)(m);

		}
	}

	coefs.save("coefs_arma.txt",arma::raw_ascii);



	return 0;
}
