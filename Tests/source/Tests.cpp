

#include "Tests.hpp"


#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>
#include <Constants.hpp>
#include <RigidBodyKinematics.hpp>

#include <vtkObjectFactory.h>
#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATMassProperties.hpp>

#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkPolyDataNormals.h>
#include <RigidBodyKinematics.hpp>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkCellData.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <vtkPolyData.h>
#include <assert.h>
#include <vtkTriangleFilter.h>
#include <vtkCleanPolyData.h>

void TestsSBCore::run() {


	TestsSBCore::test_sbgat_pgm();
	TestsSBCore::test_sbgat_mass_properties();

	TestsSBCore::test_loading_shape();
	TestsSBCore::test_geometrical_measures();
	TestsSBCore::test_pgm_consistency_cube();
	TestsSBCore::test_pgm_consistency_ellipsoid();
	TestsSBCore::test_spherical_harmonics_consistency();
	TestsSBCore::test_spherical_harmonics_invariance();


}

void TestsSBCore::test_sbgat_mass_properties(){



}


/**
This check ensures that VTK properly computes the facet normals of a cube
*/
void TestsSBCore::test_sbgat_pgm() {

	std::cout << "- Running test_sbgat_pgm ..." << std::endl;


	{
		vtkSmartPointer<vtkSphereSource> source = 
		vtkSmartPointer<vtkSphereSource>::New();

		source->SetCenter(0.0, 0.0, 0.0);
		source->SetPhiResolution(100);
		source->SetThetaResolution(100);

		vtkSmartPointer<vtkTriangleFilter> triangleFilter =
		vtkSmartPointer<vtkTriangleFilter>::New();
		triangleFilter -> SetInputConnection(source->GetOutputPort());
		triangleFilter -> Update();

		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = 
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputConnection(triangleFilter->GetOutputPort());
		cleanPolyData->Update();

		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();

		pgm_filter -> SetInputConnection(cleanPolyData -> GetOutputPort());
		pgm_filter -> Update();




	}



	{
		vtkSmartPointer<vtkCubeSource> source = 
		vtkSmartPointer<vtkCubeSource>::New();

		source->SetCenter(0.0, 0.0, 0.0);


		vtkSmartPointer<vtkTriangleFilter> triangleFilter =
		vtkSmartPointer<vtkTriangleFilter>::New();
		triangleFilter -> SetInputConnection(source->GetOutputPort());
		triangleFilter -> Update();

		vtkSmartPointer<vtkCleanPolyData> cleanPolyData = 
		vtkSmartPointer<vtkCleanPolyData>::New();
		cleanPolyData->SetInputConnection(triangleFilter->GetOutputPort());
		cleanPolyData->Update();


		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter -> SetInputConnection(cleanPolyData -> GetOutputPort());
		pgm_filter -> Update();

		vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();
		mass_filter -> SetInputConnection(cleanPolyData -> GetOutputPort());
		mass_filter -> Update();
		

		double * bounds =  mass_filter -> GetBoundingBox();

		std::cout << bounds[0] << " "  << bounds[1] << " " << bounds[2] << " " << bounds[3] << " " << bounds[4] << " " << bounds[5] << "\n";

		// Points that should be outside
		double p0[3] = {0.7,0,0};
		double p1[3] = {0.50000001,0,0};
		
		// Points that should be inside
		double p2[3] = {0.4,0,0};
		double p3[3] = {0.49999999,0,0};
		
		assert(!pgm_filter -> Contains(p0));
		assert(!pgm_filter -> Contains(p1));
		assert(pgm_filter -> Contains(p2));
		assert(pgm_filter -> Contains(p3));


	}



	std::cout << "- Done running test_sbgat_pgm ..." << std::endl;







}


/**
This check ensures that the loaded shape
model has the proper number of facets, vertices and edges
*/
void TestsSBCore::test_loading_shape() {

	std::cout << "- Running test_loading_shape ..." << std::endl;

	SBGAT_CORE::ShapeModel shape_model("", nullptr);
	SBGAT_CORE::ShapeModelImporter shape_io("../cube.obj", 1);
	shape_io.load_shape_model(&shape_model,false);

	assert(shape_model.get_NFacets() == 12);
	assert(shape_model.get_NVertices() == 8);
	assert(shape_model.get_NEdges() == 18);

	std::cout << "-- test_loading_shape successful" << std::endl;

}




/**
This check ensures that the geometrical measures such as volume 
 surface area, center of mass and inertia moments are properly computed
*/
void TestsSBCore::test_geometrical_measures() {

	std::cout << "- Running test_geometrical_measures ..." << std::endl;

	SBGAT_CORE::ShapeModel shape_model("", nullptr);
	SBGAT_CORE::ShapeModelImporter shape_io("../cube.obj", 0.5);
	shape_io.load_shape_model(&shape_model,false);


	arma::vec cube_moments = {1,1,1};
	arma::vec cube_com = {1,1,1};

	cube_moments /= 6.;

	cube_com /= 4.;

	shape_model . update_mass_properties();

	arma::mat inertia = shape_model . get_inertia();

	assert(std::abs(1./8 - shape_model.get_volume()) < 1e-10);
	assert(std::abs(6./4 - shape_model.get_surface_area()) < 1e-10);
	assert(arma::norm(cube_com - shape_model.get_center_of_mass()) < 1e-10);
	assert(arma::norm(cube_moments - arma::eig_sym(inertia)) / arma::norm(cube_moments) < 1e-10);


	std::cout << "-- test_geometrical_measures successful" << std::endl;

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
		shape_io.load_shape_model(&shape_model,true);
		SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
		double mu = shape_model.get_volume() * 1e6 * arma::datum::G;
		arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

		assert(arma::norm(acc_true - acc) / arma::norm(acc) < 1e-12);
	}

	{
		SBGAT_CORE::ShapeModel shape_model("", nullptr);
		SBGAT_CORE::ShapeModelImporter shape_io("../cube_50k.obj", 1);
		shape_io.load_shape_model(&shape_model,true);
		double mu = shape_model.get_volume() * 1e6 * arma::datum::G;




		SBGAT_CORE::DynamicAnalyses dyn_an(&shape_model);
		arma::vec acc = dyn_an.pgm_acceleration(X.colptr(0), mu);

		assert(arma::norm(acc_true - acc) / arma::norm(acc) < 1e-12);
	}

	{
		SBGAT_CORE::ShapeModel shape_model("", nullptr);
		SBGAT_CORE::ShapeModelImporter shape_io("../cube_200k.obj", 1);
		shape_io.load_shape_model(&shape_model,true);
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
	shape_io.load_shape_model(&shape_model,false);

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

		shape_io.load_shape_model(&eros,false);
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

		shape_io.load_shape_model(&eros,false);
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









