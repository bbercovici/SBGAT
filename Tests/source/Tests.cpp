

#include "Tests.hpp"

void TestsSBCore::run() {
	TestsSBCore::test_loading_shape();
	TestsSBCore::test_pgm_consistency_cube();
	TestsSBCore::test_pgm_consistency_ellipsoid();
	TestsSBCore::test_spherical_harmonics_consistency();
	TestsSBCore::test_spherical_harmonics_invariance();



}


/**
This check ensures that the loaded shape
model has the proper number of facets, vertices and edges
*/
void TestsSBCore::test_loading_shape() {

	std::cout << "- Running test_loading_shape ..." << std::endl;

	SBGAT_CORE::ShapeModel shape_model("", nullptr);
	SBGAT_CORE::ShapeModelImporter shape_io("../cube.obj", 1);
	shape_io.load_shape_model(&shape_model);

	assert(shape_model.get_NFacets() == 12);
	assert(shape_model.get_NVertices() == 8);
	assert(shape_model.get_NEdges() == 18);

	std::cout << "-- test_loading_shape successful" << std::endl;

}


/**
This test checks that the Polyhedron Gravity Model
computed around a cubic shape is not depending upon
the cube resolution and consistent with the analytical expression
*/
void TestsSBCore::test_pgm_consistency_cube() {

	std::cout << "- Running test_pgm_consistency_cube ..." << std::endl;


	// The attracting shape is a cube of dimensions 1 x 1 x 1 m of density rhp = 1e6 kg/m^3
	// The analytic acceleration at (1,2,3) (m) in the shape model's barycentric frame is computed
	// Assumes that G = 6.67408e-11 m^3 / (kg * s ^2)
	arma::vec X = {1, 2, 3};
	arma::vec acc_true = {
		-1.273782722739791e-06,
		-2.548008881415967e-06,
		-3.823026510474731e-06
	};


	// Shape models of increasing resolutions but still representative
	// of a cube are loaded. They should all yield the same acceleration
	// at the specified query point

	{
		SBGAT_CORE::ShapeModel shape_model("", nullptr);
		SBGAT_CORE::ShapeModelImporter shape_io("../cube.obj", 1);
		shape_io.load_shape_model(&shape_model);
		SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
		double mu = shape_model.get_volume() * 1e6 * arma::datum::G;
		arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

		assert(arma::norm(acc_true - acc) / arma::norm(acc) < 1e-12);
	}

	{
		SBGAT_CORE::ShapeModel shape_model("", nullptr);
		SBGAT_CORE::ShapeModelImporter shape_io("../cube_50k.obj", 1);
		shape_io.load_shape_model(&shape_model);
		double mu = shape_model.get_volume() * 1e6 * arma::datum::G;




		SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
		arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

		assert(arma::norm(acc_true - acc) / arma::norm(acc) < 1e-12);
	}

	{
		SBGAT_CORE::ShapeModel shape_model("", nullptr);
		SBGAT_CORE::ShapeModelImporter shape_io("../cube_200k.obj", 1);
		shape_io.load_shape_model(&shape_model);
		double mu = shape_model.get_volume() * 1e6 * arma::datum::G;


		SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
		arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

		assert(arma::norm(acc_true - acc) / arma::norm(acc) < 1e-12);
	}

	std::cout << "-- test_pgm_consistency_cube successful" << std::endl;


}



/**
This test ensures that the Polyhedron Gravity Model
computed around an ellispoid shape is consistent
with the analytical expression
*/
void TestsSBCore::test_pgm_consistency_ellipsoid() {

	std::cout << "- Running test_pgm_consistency_ellipsoid ..." << std::endl;

	// The attracting shape is an ellipsoid of semi-major axes 3 x 2 x 1 m of density rho = 1e6 kg/m^3
	// The analytic acceleration at (1,2,3) (m) in the shape model's barycentric frame is computed
	// Assumes that G = 6.67408e-11 m^3 / (kg * s ^2)
	arma::vec X = {1, 2, 3};
	arma::vec acc_true = {
		-2.19160852e-05,  
		-5.18364044e-05,  
		-8.79434337e-05
	};

	SBGAT_CORE::ShapeModel shape_model("", nullptr);
	SBGAT_CORE::ShapeModelImporter shape_io("../ellipsoid.obj", 1);
	shape_io.load_shape_model(&shape_model);

	SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
	double mu = shape_model.get_volume() * 1e6 * arma::datum::G;

	arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

	assert(arma::norm(acc_true - acc) / arma::norm(acc) < 5e-5);

	std::cout << "-- test_pgm_consistency_ellipsoid successful" << std::endl;

}


/**
This test checks the consistency of the spherical harmonics coefficients computation 
about a shape model of Eros against

*/

void TestsSBCore::test_spherical_harmonics_consistency() {

	std::cout << "- Running test_spherical_harmonics_consistency ..." << std::endl;


	SBGAT_CORE::FrameGraph frame_graph;
	SBGAT_CORE::ShapeModel eros("B", &frame_graph);
	SBGAT_CORE::ShapeModelImporter shape_io("../eros_64.obj", 1);

	shape_io.load_shape_model(&eros,false);
	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&eros);

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of Eros (kg/km^3). The input shape model has 
	// its coordinates expressed in km 
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

	// The standard gravitational parameter of Eros is computed
	double mu = eros.get_volume() * (density * arma::datum::G / 1e9);

	// The accelerations are evaluated at the query point
	arma::vec pos = {10,10,10};

	// Accelerations obtained from the legacy MEX implementation 
	// of this spherical harmonics code by Yu Takahashi
	// and Siamak Hesar

	arma::vec acc_true = { -0.567824977098112e-6,
		-0.846742135468519e-6,
		-0.836664679681573e-6
	};

	arma::vec acc_sph = dynamic_analyses.spherical_harmo_acc(degree,
		ref_radius,
		mu,
		pos, 
		Cnm_total,
		Snm_total);

	assert(arma::norm(acc_true - acc_sph) / arma::norm(acc_sph) < 1e-10);

	std::cout << "-- test_spherical_harmonics_consistency successful" << std::endl;

}

void TestsSBCore::test_spherical_harmonics_invariance(){

	std::cout << "- Running test_spherical_harmonics_invariance ..." << std::endl;

	arma::vec acc_sph_km(3);
	arma::vec acc_sph_m(3);


	{

		SBGAT_CORE::FrameGraph frame_graph;
		SBGAT_CORE::ShapeModel eros("B", &frame_graph);
		SBGAT_CORE::ShapeModelImporter shape_io("../eros_64.obj", 1);

		shape_io.load_shape_model(&eros);
		SBGAT_CORE::DynamicAnalyses dynamic_analyses(&eros);

		int degree = 10;

		double density = 2670000000000.0;

		double ref_radius = 16;

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

		arma::vec pos = {10,10,10};

		double mu = eros.get_volume() * (density * arma::datum::G / 1e9);

		acc_sph_km = dynamic_analyses.spherical_harmo_acc(degree,
			ref_radius,
			mu,
			pos, 
			Cnm_total,
			Snm_total);
	}

	{

		SBGAT_CORE::FrameGraph frame_graph;
		SBGAT_CORE::ShapeModel eros("B", &frame_graph);
		SBGAT_CORE::ShapeModelImporter shape_io("../eros_64.obj", 1000);

		shape_io.load_shape_model(&eros);
		SBGAT_CORE::DynamicAnalyses dynamic_analyses(&eros);

		int degree = 10;

		double density = 2670;

		double ref_radius = 16000;

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

		arma::vec pos = {10000,10000,10000};

		double mu = eros.get_volume() * (density * arma::datum::G);

		acc_sph_m = dynamic_analyses.spherical_harmo_acc(degree,
			ref_radius,
			mu,
			pos, 
			Cnm_total,
			Snm_total);
	}

	assert(arma::norm(acc_sph_m - 1e3 * acc_sph_km)/arma::norm(acc_sph_km) < 1e-10);


	std::cout << "-- test_spherical_harmonics_invariance successful" << std::endl;


}









