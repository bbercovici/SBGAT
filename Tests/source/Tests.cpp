

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
#include <vtkOBJReader.h>
#include <vtkCellCenters.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkLinearSubdivisionFilter.h>

void TestsSBCore::run() {

	TestsSBCore::test_sbgat_mass_properties();
	TestsSBCore::test_sbgat_pgm();
	TestsSBCore::test_sbgat_pgm_speed();

	// TODO: reimplement in SbgatCore
	// TestsSBCore::test_spherical_harmonics_consistency();
	// TestsSBCore::test_spherical_harmonics_invariance();


}


/**
This test compares some geometric measures as computed by SBGAT 
to analytical values
*/
void TestsSBCore::test_sbgat_mass_properties(){


	vtkSmartPointer<vtkCubeSource> source = 
	vtkSmartPointer<vtkCubeSource>::New();
	source -> SetCenter(0.0, 0.0, 0.0);

	vtkSmartPointer<vtkTriangleFilter> triangleFilter =
	vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter -> SetInputConnection(source->GetOutputPort());
	triangleFilter -> Update();

	vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();

	subdivisionFilter -> SetInputConnection(triangleFilter -> GetOutputPort());
	subdivisionFilter -> SetNumberOfSubdivisions(6);
	subdivisionFilter -> Update();

	vtkSmartPointer<vtkCleanPolyData> cleanPolyData = 
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPolyData->SetInputConnection(subdivisionFilter->GetOutputPort());
	cleanPolyData->Update();

	vtkSmartPointer<vtkTransform> transform_rot = vtkSmartPointer<vtkTransform>::New();
	transform_rot->RotateWXYZ(10, 0, 1, 0);
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter_rot = 
	vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter_rot->SetTransform(transform_rot);
	transformFilter_rot->SetInputConnection(cleanPolyData->GetOutputPort());
	transformFilter_rot->Update();


	vtkSmartPointer<vtkTransform> transform_trans = vtkSmartPointer<vtkTransform>::New();
	arma::vec translation_vector = { 1, 2, 3};
	transform_trans->Translate (translation_vector.colptr(0));
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter_trans = 
	vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter_trans->SetTransform(transform_trans);
	transformFilter_trans->SetInputConnection(transformFilter_rot->GetOutputPort());
	transformFilter_trans->Update();

	vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();
	mass_filter -> SetInputConnection(transformFilter_trans -> GetOutputPort());
	mass_filter -> Update();

	// the inertia moments are invariant by rotation/translation
	arma::vec inertia_moments = {1./6,1./6,1./6}; 
	auto inertia_moments_sbgat = mass_filter -> GetInertiaMoments();

	// The center of mass should be right where the translation put it 
	auto com_sbgat = mass_filter -> GetCenterOfMass();
	
	assert(arma::norm(inertia_moments - inertia_moments_sbgat)/arma::norm(inertia_moments_sbgat) < 1e-8);
	assert(arma::norm(translation_vector - com_sbgat)/arma::norm(com_sbgat) < 1e-8);

}


/**
This test computes the surface accelerations at the center of each facet over a polydata
of Eros for benchmarking purposes
*/
void TestsSBCore::test_sbgat_pgm_speed(){
	std::cout << "- Running test_sbgat_pgm_speed ..." << std::endl;
	
	// Reading
	vtkNew<vtkOBJReader> reader;
	reader -> SetFileName("../eros_64.obj");
	reader -> Update(); 


	vtkSmartPointer<vtkTriangleFilter> triangleFilter =
	vtkSmartPointer<vtkTriangleFilter>::New();
	triangleFilter -> SetInputConnection(reader->GetOutputPort());
	triangleFilter -> Update();

	vtkSmartPointer<vtkCleanPolyData> cleanPolyData = 
	vtkSmartPointer<vtkCleanPolyData>::New();
	cleanPolyData->SetInputConnection(triangleFilter->GetOutputPort());
	cleanPolyData->Update();

	std::cout << "-- Creating dyads...\n";
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleanPolyData -> GetOutputPort());
	pgm_filter -> Update();
	std::cout << "-- Done creating dyads...\n";


	vtkSmartPointer<vtkPolyData> polydata = cleanPolyData -> GetOutput();

	vtkSmartPointer<vtkCellCenters> cellCentersFilter = 
	vtkSmartPointer<vtkCellCenters>::New();
	
	cellCentersFilter -> SetInputConnection(cleanPolyData -> GetOutputPort());
	cellCentersFilter -> Update();

	arma::mat surface_accelerations(cellCentersFilter -> GetOutput() -> GetNumberOfPoints(),3);

	assert(polydata -> GetNumberOfCells() == cellCentersFilter -> GetOutput() -> GetNumberOfPoints());
	auto start = std::chrono::system_clock::now();
	std::cout << "-- Computing pgm accelerations at " << cellCentersFilter -> GetOutput() -> GetNumberOfPoints() << " facet centers over the surface of Eros. This may take a few minutes ...\n";
	for (vtkIdType i = 0; i < cellCentersFilter -> GetOutput() -> GetNumberOfPoints(); ++i){
		double p[3];
		cellCentersFilter -> GetOutput() -> GetPoint(i, p);
		surface_accelerations.row(i) = pgm_filter -> ComputePgmAcceleration(p ,2670.).t();
	}
	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "-- Done computing pgm accelerations in " << elapsed_seconds.count() << " s\n";

	std::cout << "- Done running test_sbgat_pgm_speed ..." << std::endl;

}


/**
This test computes the pgm acceleration and potential
about a cube and compares the SBGAT outputs to analytical ones 
*/
void TestsSBCore::test_sbgat_pgm() {

	std::cout << "- Running test_sbgat_pgm ..." << std::endl;

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



	// The attracting shape is a cube of dimensions 1 x 1 x 1 m of density rho = 1e6 kg/m^3
	// The analytic potential and acceleration at (1,2,3) (m) in the shape model's barycentric frame is computed
	// Assumes that G = 6.67408e-11 m^3 / (kg * s ^2)


		// Queried point for pgm validation
		double p4[3] = {1, 2, 3};

		arma::vec acc_true = {
			-1.273782722739791e-06,
			-2.548008881415967e-06,
			-3.823026510474731e-06
		};
		double pot_true = 0.26726619638669064 * arma::datum::G * 1e6;
		arma::vec pgm_acc = pgm_filter -> ComputePgmAcceleration(p4,1e6);
		double pgm_pot = pgm_filter -> ComputePgmPotential(p4,1e6);

		assert(arma::norm(pgm_acc - acc_true)/arma::norm(acc_true) < 1e-10);
		assert(std::abs(pgm_pot - pot_true)/std::abs(pot_true) < 1e-10);



	}



	std::cout << "- Done running test_sbgat_pgm ..." << std::endl;




}


/**
This test checks the consistency of the spherical harmonics coefficients computation 
about a shape model of Eros against, evaluates the resulting acceleration, and compares
it to one extracted from the polyhedron gravity model acceleration

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









