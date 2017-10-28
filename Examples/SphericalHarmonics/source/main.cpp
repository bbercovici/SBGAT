/** 
@file   main.cpp
@Author Benjamin Bercovici (bebe0705@colorado.edu)
@date   September, 2017
@brief  Example demonstrating the computation of gravity spherical harmonics
assuming a constant density from an input polyhedral shape model
*/

#include <iostream>
#include <armadillo>
#include <chrono>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>

int main( int argc, char** argv ) {


	SBGAT_CORE::FrameGraph frame_graph;
	SBGAT_CORE::ShapeModel eros("B", &frame_graph);
	SBGAT_CORE::ShapeModelImporter shape_io("../eros_64.obj", 1);

	shape_io.load_shape_model(&eros);
	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&eros);

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of Eros (kg/km^3)
	double density = 2670000000000.0;

	// Reference radius of Eros (km)
	double ref_radius = 16;

	// Flag set to true for normalized coefficients, false for unnormalized ones
	bool normalized = true;

	arma::mat Cnm_total;
	arma::mat Snm_total;

	dynamic_analyses.compute_exterior_sh_coefs(
		Cnm_total,
		Snm_total,
		degree,
		ref_radius,
		density,
		normalized);

	// This coefficients can then be used to evaluate an acceleration at a query point
	// outside of the sphere circumscribing the input shape
	arma::vec pos = {10,10,10};

	// The standard gravitational parameter of Eros is computed
	double mu = eros.get_volume() * (density * arma::datum::G / 1e9);

	arma::vec acc_sph = dynamic_analyses.spherical_harmo_acc(degree,
		ref_radius,
		mu,
		pos, 
		Cnm_total,
		Snm_total);

	std::cout << "Spherical harmonics acceleration at the query point (km/s^2): " << std::endl;
	std::cout << acc_sph << std::endl;

	arma::vec acc_pgm = dynamic_analyses.pgm_acceleration(pos.colptr(0) , density);

	std::cout << "Polyhedron gravity model acceleration at the query point (km/s^2): " << std::endl;
	std::cout << acc_pgm << std::endl;


	// The coefficients are stored in a more convenient form 
	arma::mat coefs = arma::zeros<arma::mat>((degree + 1) * (degree + 2)/2 - 1,4);
	
	for (int n = 1 ; n < degree + 1 ; ++n){
		for (int m = 0 ; m <= n; ++m){

			int index = (n) * (n + 1 )/2 + m - 1;

			coefs.row(index)(0) = n;
			coefs.row(index)(1) = m;
			coefs.row(index)(2) = Cnm_total(n,m);
			coefs.row(index)(3) = Snm_total(n,m);

		}
	}

	std::cout << "Cbar coefficients : " << std::endl;
	std::cout << Cnm_total << std::endl;


	std::cout << "Sbar coefficients : " << std::endl;
	std::cout << Snm_total << std::endl;


	// coefs now holds the normalized spherical harmonics coefficients
	// the matrix's 4 columns are labelled as follows:
	// degree n -- order m -- Cnm -- Snm

	std::cout << "Sorted coefficients (degree n -- order m -- Cnm -- Snm) : "<< std::endl;
	std::cout << coefs << std::endl;

	// The sorted coefficient table is saved to a file
	coefs.save("eros_spherical_coords_normalized.txt",arma::raw_ascii);

	// The unsorted coefficients are sorted as well	
	Cnm_total.save("eros_Cbar.txt",arma::raw_ascii);
	Snm_total.save("eros_Sbar.txt",arma::raw_ascii);




	return 0;
}
