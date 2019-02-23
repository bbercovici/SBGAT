/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "Tests.hpp"

#include <vtkObjectFactory.h>
#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATMassProperties.hpp>
#include <SBGATSphericalHarmo.hpp>

#include <SBGATObsRadar.hpp>
#include <SBGATObsLightcurve.hpp>
#include <SBGATFrameGraph.hpp>
#include <SBGATShapeUncertainty.hpp>
#include <SBGATTransformShape.hpp>
#include <SBGATObjWriter.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>

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
#include <boost/progress.hpp>


void TestsSBCore::run() {	
	TestsSBCore::test_sbgat_transform_shape();
	TestsSBCore::test_frame_conversion();
	TestsSBCore::test_sbgat_mass_properties();
	TestsSBCore::test_sbgat_pgm_cube();
	TestsSBCore::test_sbgat_pgm_sphere();
	TestsSBCore::test_sbgat_pgm_speed();
	TestsSBCore::test_spherical_harmonics_coefs_consistency();
	TestsSBCore::test_spherical_harmonics_partials_consistency();
	TestsSBCore::test_sbgat_shape_uq();
	TestsSBCore::test_PGM_UQ_partials();
	TestsSBCore::test_PGM_UQ_covariance_consistency();
	TestsSBCore::test_PGM_UQ_cube();
	TestsSBCore::test_PGM_UQ_itokawa_m();
	TestsSBCore::test_PGM_UQ_itokawa_km();


	// TestsSBCore::test_lightcurve_obs();
	// TestsSBCore::test_radar_obs();


	std::cout << "All tests passed.\n";

}




/**
This test compares some geometric measures as computed by SBGAT 
to analytical values
*/
void TestsSBCore::test_sbgat_mass_properties(){

	std::cout << "- Running test_sbgat_mass_properties ..." << std::endl;

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
	auto com_sbgat = mass_filter -> GetCenterOfMass();

	assert(mass_filter -> CheckClosed());

	// the inertia moments are invariant by rotation/translation
	arma::vec inertia_moments = {1./6,1./6,1./6}; 
	double r_avg = std::cbrt(3./ (4 * arma::datum::pi));
	inertia_moments = inertia_moments / (r_avg * r_avg) + arma::eig_sym(RBK::tilde(translation_vector) * RBK::tilde(translation_vector).t())/ (r_avg * r_avg);
	auto inertia_moments_sbgat = mass_filter -> GetNormalizedInertiaMoments();

	assert(arma::norm(inertia_moments - inertia_moments_sbgat)/arma::norm(inertia_moments_sbgat) < 1e-8);
	assert(arma::norm(translation_vector - com_sbgat)/arma::norm(com_sbgat) < 1e-8);


	double volume_itokawa_m,volume_itokawa_km;



	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName("../../resources/shape_models/itokawa_8.obj");
	reader -> Update(); 
	mass_filter -> SetInputConnection(reader -> GetOutputPort());
	mass_filter -> SetScaleKiloMeters();
	mass_filter -> Modified();
	mass_filter -> Update();

	volume_itokawa_km = mass_filter -> GetVolume();

	reader -> SetFileName("../../resources/shape_models/itokawa_8_scaled.obj");
	reader -> Update(); 
	mass_filter -> SetInputConnection(reader -> GetOutputPort());
	mass_filter -> SetScaleMeters();
	mass_filter -> Modified();
	mass_filter -> Update();
	volume_itokawa_m = mass_filter -> GetVolume();

	assert(std::abs(volume_itokawa_m - volume_itokawa_km) / volume_itokawa_m < 1e-7);

	std::cout << "- Done running test_sbgat_mass_properties" << std::endl;

}


/**
This test computes the surface accelerations at the center of each facet over a polydata
of KW4 for benchmarking purposes
*/
void TestsSBCore::test_sbgat_pgm_speed(){
	std::cout << "- Running test_sbgat_pgm_speed ..." << std::endl;
	
	std::string filename  = "../../resources/shape_models/KW4Alpha.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Creating the PGM dyads
	std::cout << "-- Creating dyads...\n";
	double density = 2000;

	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(density); // density in kg/m^3
	pgm_filter -> SetScaleKiloMeters();
	pgm_filter -> Update();
	std::cout << "-- Done creating dyads...\n";

	// Creating a filter that will extract the center of each facet
	vtkSmartPointer<vtkPolyData> polydata = reader -> GetOutput();
	vtkSmartPointer<vtkCellCenters> cellCentersFilter = 
	vtkSmartPointer<vtkCellCenters>::New();
	cellCentersFilter -> SetInputConnection(reader -> GetOutputPort());
	cellCentersFilter -> Update();

	// Preallocating
	arma::mat surface_accelerations(3,cellCentersFilter -> GetOutput() -> GetNumberOfPoints());

	assert(polydata -> GetNumberOfCells() == cellCentersFilter -> GetOutput() -> GetNumberOfPoints());

	int N_f = polydata -> GetNumberOfCells();
	std::cout << "-- Computing pgm accelerations at " << cellCentersFilter -> GetOutput() -> GetNumberOfPoints() << " facet centers over the surface of " << filename <<  " . This may take a few minutes ...\n";
	
	boost::progress_display progress(N_f);
	auto start = std::chrono::system_clock::now();
	#pragma omp parallel for
	for (vtkIdType i = 0; i < N_f; ++i){

		// Input position must be in meters

		arma::vec::fixed<3> p = 1000 * pgm_filter -> GetFacetCenter(i);

		surface_accelerations.col(i) = pgm_filter -> GetAcceleration(p);
		++progress;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "\n-- Done computing pgm accelerations in " << elapsed_seconds.count() << " s\n";

	std::vector<double> slopes;
	std::vector<double> inertial_potentials;
	std::vector<double> body_fixed_potentials;
	std::vector<double> inertial_acc_magnitudes;
	std::vector<double> body_fixed_acc_magnitudes;
	arma::vec::fixed<3> omega = {0,0,2 * arma::datum::pi / (2.7650 * 3600)};


	std::vector<unsigned int> queried_elements;
	for (int i = 0; i < polydata -> GetNumberOfCells() ; ++i){
		queried_elements.push_back(i);
	}

	start = std::chrono::system_clock::now();

	SBGATPolyhedronGravityModel::ComputeSurfacePGM(polydata,queried_elements,false,density,omega,
		slopes,
		inertial_potentials,
		body_fixed_potentials,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes);

	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	std::cout << "-- Done computing pgm slopes, potentials and accelerations in " << elapsed_seconds.count() << " s\n";

	assert(N_f == inertial_acc_magnitudes.size());

	for (int i = 0; i < N_f; ++i){
		assert(std::abs(arma::norm(surface_accelerations.col(i)) - inertial_acc_magnitudes[i]) / inertial_acc_magnitudes[i] < 1e-5);
	}

	
	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();

	mass_properties -> SetInputData(polydata);
	mass_properties -> SetScaleKiloMeters();

	mass_properties -> Update();
	double mass = mass_properties -> GetVolume() * density;

	SBGATPolyhedronGravityModel::SaveSurfacePGM(polydata,
		queried_elements,
		false,
		mass,
		omega,
		slopes,
		inertial_potentials,
		body_fixed_potentials,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes,
		"../output/kw4_pgm.json");



	std::cout << "- Done running test_sbgat_pgm_speed" << std::endl;

}


/**
This test computes the pgm acceleration and potential
about a unit cube and compares the SBGAT outputs to analytical ones 
*/
void TestsSBCore::test_sbgat_pgm_cube() {

	std::cout << "- Running test_sbgat_pgm_cube ..." << std::endl;

	{
		vtkSmartPointer<vtkCubeSource> source = 
		vtkSmartPointer<vtkCubeSource>::New();

		source -> SetCenter(0.0, 0.0, 0.0);

		double density = 1e6;

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


		pgm_filter -> SetDensity(density);
		pgm_filter -> SetScaleMeters();
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
		arma::vec p4 = {1, 2, 3};

		arma::vec acc_true = {
			-1.273782722739791e-06,
			-2.548008881415967e-06,
			-3.823026510474731e-06
		};
		double pot_true = 0.26726619638669064 * arma::datum::G * density;


		arma::vec pgm_acc = pgm_filter -> GetAcceleration(p4.colptr(0));
		double pgm_pot = pgm_filter -> GetPotential(p4.colptr(0));

		assert(arma::norm(pgm_acc - acc_true)/arma::norm(acc_true) < 1e-10);
		assert(std::abs(pgm_pot - pot_true)/std::abs(pot_true) < 1e-10);

		pgm_acc = pgm_filter -> GetAcceleration(p4);
		pgm_pot = pgm_filter -> GetPotential(p4);

		double pgm_pot2;
		arma::vec::fixed<3> pgm_acc2;
		arma::mat::fixed<3,3> pgm_gravity_gradient_matrix;
		pgm_filter -> GetPotentialAccelerationGravityGradient(p4,pgm_pot2, pgm_acc2,pgm_gravity_gradient_matrix);

		assert(arma::norm(pgm_acc - acc_true)/arma::norm(acc_true) < 1e-10);
		assert(std::abs(pgm_pot - pot_true)/std::abs(pot_true) < 1e-10);

		assert(arma::norm(pgm_acc2 - acc_true)/arma::norm(acc_true) < 1e-10);
		assert(std::abs(pgm_pot2 - pot_true)/std::abs(pot_true) < 1e-10);

		arma::vec::fixed<3> p4_perturbed = (1 + 1e-3) * p4;

		arma::vec::fixed<3> pgm_acc2_perturbed = pgm_filter -> GetAcceleration(p4_perturbed);

		auto d_acc_true = pgm_acc2_perturbed - pgm_acc2;

		arma::vec::fixed<3> d_acc_linear = pgm_gravity_gradient_matrix * (p4_perturbed - p4);

		assert(100 * arma::norm(d_acc_linear - d_acc_true)/arma::norm(d_acc_true) < 1e0);

	}

	std::cout << "- Done running test_sbgat_pgm_cube" << std::endl;

}


/**
This test computes the pgm acceleration and potential
about a sphere and compares the SBGAT outputs to analytical ones 
*/
void TestsSBCore::test_sbgat_pgm_sphere() {

	std::cout << "- Running test_sbgat_pgm_sphere ..." << std::endl;

	
	vtkSmartPointer<vtkSphereSource> source = 
	vtkSmartPointer<vtkSphereSource>::New();

	source -> SetCenter(0.0, 0.0, 0.0);
	source -> SetRadius(6370);
	source -> SetThetaResolution (300);
	source -> SetPhiResolution (300);

	double density = 5.51e3;

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
	pgm_filter -> SetDensity(density);
	pgm_filter -> SetScaleKiloMeters();
	pgm_filter -> Update();

	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();

	mass_properties -> SetInputConnection(cleanPolyData -> GetOutputPort());
	mass_properties -> SetScaleKiloMeters();
	mass_properties -> Update();

		// Queried point for pgm validation
	arma::vec p = {6371e3, 0, 0};

	arma::vec acc_true = - arma::datum::G * density * mass_properties -> GetVolume() / arma::dot(p,p) * arma::normalise(p);
	double pot_true = arma::datum::G * density * mass_properties -> GetVolume() / arma::norm(p);

	arma::vec pgm_acc = pgm_filter -> GetAcceleration(p);
	double pgm_pot = pgm_filter -> GetPotential(p);


	assert(arma::norm(pgm_acc - acc_true)/arma::norm(acc_true) < 1e-4);
	assert(std::abs(pgm_pot - pot_true)/std::abs(pot_true) < 1e-4);

	pgm_acc = pgm_filter -> GetAcceleration(p);
	pgm_pot = pgm_filter -> GetPotential(p);

	double pgm_pot2;
	arma::vec::fixed<3> pgm_acc2;
	arma::mat::fixed<3,3> pgm_gravity_gradient_matrix;
	pgm_filter -> GetPotentialAccelerationGravityGradient(p,pgm_pot2, pgm_acc2,pgm_gravity_gradient_matrix);

	assert(arma::norm(pgm_acc - acc_true)/arma::norm(acc_true) < 1e-4);
	assert(std::abs(pgm_pot - pot_true)/std::abs(pot_true) < 1e-4);

	assert(arma::norm(pgm_acc2 - acc_true)/arma::norm(acc_true) < 1e-4);
	assert(std::abs(pgm_pot2 - pot_true)/std::abs(pot_true) < 1e-4);

	arma::vec::fixed<3> p4_perturbed = (1 + 1e-3) * p;

	arma::vec::fixed<3> pgm_acc2_perturbed = pgm_filter -> GetAcceleration(p4_perturbed);

	auto d_acc_true = pgm_acc2_perturbed - pgm_acc2;

	arma::vec::fixed<3> d_acc_linear = pgm_gravity_gradient_matrix * (p4_perturbed - p);

	assert(100 * arma::norm(d_acc_linear - d_acc_true)/arma::norm(d_acc_true) < 1e0);

	

	std::cout << "- Done running test_sbgat_pgm_sphere" << std::endl;

}


/**
This test checks the consistency of the spherical harmonics coefficients computation 
about a shape model of KW4 
*/

void TestsSBCore::test_spherical_harmonics_coefs_consistency() {

	std::cout << "- Running test_spherical_harmonics_coefs_consistency ..." << std::endl;


	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName("../../resources/shape_models/KW4Alpha.obj");
	reader -> Update(); 

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of KW4 (kg/m^3). 
	double density = 2000.0;

	// Reference radius of KW4 (m)
	double ref_radius = 1.317/2 * 1000;

	// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(density);
	pgm_filter -> SetScaleKiloMeters();// The input shape model has its coordinates expressed in km 
	pgm_filter -> Update();

	// An instance of SBGATSphericalHarmo is created to compute and evaluate the spherical 
	// expansion of the gravity field about the considered shape model
	vtkSmartPointer<SBGATSphericalHarmo> spherical_harmonics = vtkSmartPointer<SBGATSphericalHarmo>::New();
	spherical_harmonics -> SetInputConnection(reader -> GetOutputPort());
	spherical_harmonics -> SetDensity(density);
	spherical_harmonics -> SetScaleKiloMeters();
	spherical_harmonics -> SetReferenceRadius(ref_radius);
	spherical_harmonics -> IsNormalized(); // can be skipped as normalized coefficients is the default parameter
	spherical_harmonics -> SetDegree(degree);
	spherical_harmonics -> Update();


	// The spherical harmonics are saved to a file
	spherical_harmonics -> SaveToJson("../output/KW4Alpha_harmo.json");

	
	// The accelerations are evaluated at the query point
	arma::vec::fixed<3> pos = {3e3,5e3,-2e3};
	arma::vec::fixed<3> pgm_acc = pgm_filter -> GetAcceleration(pos);
	arma::vec::fixed<3> sharm_acc = spherical_harmonics -> GetAcceleration(pos);

	// The spherical harmonics are read from the just-saved JSON file and re-evaluated
	vtkSmartPointer<SBGATSphericalHarmo> spherical_harmonics_from_file = vtkSmartPointer<SBGATSphericalHarmo>::New();
	spherical_harmonics_from_file -> LoadFromJson("../output/KW4Alpha_harmo.json");
	arma::vec::fixed<3> sharm_acc_from_file = spherical_harmonics_from_file -> GetAcceleration(pos);

	// The accelerations should be consistent with the one previously computed

	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	mass_properties -> SetInputConnection(reader -> GetOutputPort());
	mass_properties -> SetScaleKiloMeters();
	mass_properties -> Update();


	assert(arma::norm(pgm_acc - sharm_acc_from_file) / arma::norm(sharm_acc_from_file) * 100 < 1e-4);
	assert(arma::norm(sharm_acc_from_file - sharm_acc) / arma::norm(sharm_acc) * 100 < 1e-8);
	
	std::cout << "- Done running test_spherical_harmonics_coefs_consistency" << std::endl;

}



/**
This test checks the consistency of the spherical harmonics coefficients computation 
about a shape model of KW4 
*/

void TestsSBCore::test_spherical_harmonics_partials_consistency() {

	std::cout << "- Running test_spherical_harmonics_partials_consistency ..." << std::endl;

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName("../../resources/shape_models/KW4Alpha.obj");
	reader -> Update(); 

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of KW4 (kg/m^3). 
	double density = 2000.0;

	// Reference radius of KW4 (km)
	double ref_radius = 1.317/2 * 1000;

	// An instance of SBGATSphericalHarmo is created to compute and evaluate the spherical 
	// expansion of the gravity field about the considered shape model
	vtkSmartPointer<SBGATSphericalHarmo> spherical_harmonics = vtkSmartPointer<SBGATSphericalHarmo>::New();
	spherical_harmonics -> SetInputConnection(reader -> GetOutputPort());
	spherical_harmonics -> SetDensity(density);
	spherical_harmonics -> SetScaleKiloMeters();
	spherical_harmonics -> SetReferenceRadius(ref_radius);
	spherical_harmonics -> IsNormalized(); // can be skipped as normalized coefficients is the default parameter
	spherical_harmonics -> SetDegree(degree);
	spherical_harmonics -> Update();


	// The accelerations are evaluated at the query point
	arma::vec::fixed<3> pos = {3e3,5e3,-2e3};
	arma::vec::fixed<3> dpos = 1e-3 * pos;

	arma::vec::fixed<3> sharm_acc = spherical_harmonics -> GetAcceleration(pos);
	arma::vec::fixed<3> sharm_acc_pert = spherical_harmonics -> GetAcceleration(pos + dpos);
	arma::vec::fixed<3> dacc_true = sharm_acc_pert - sharm_acc;
	arma::mat::fixed<3,3> gravity_gradient_mat;
	spherical_harmonics -> GetGravityGradientMatrix(pos,gravity_gradient_mat);

	arma::vec::fixed<3> dacc_lin = gravity_gradient_mat * dpos;
	assert(arma::norm(dacc_lin - dacc_true ) / arma::norm(dacc_true) * 100 < 1e0);


	// Modifying the coefs and testing the partials
	arma::mat Cnm = spherical_harmonics-> GetCnm();
	arma::mat Snm = spherical_harmonics-> GetSnm();


	arma::mat dCnm = 1e-2 * Cnm;
	arma::mat dSnm = 1e-2 * Snm;

	dCnm(0,0) = 0;
	dSnm.col(0) *= 0;


	arma::vec dCnm_vec((degree + 1) * (degree + 2)/2);
	arma::vec dSnm_vec((degree + 1) * (degree + 2)/2 - degree - 1);

	int counter = 0;
	for (int i = 0; i <= degree; ++i){
		for (int j = 0; j <= i; ++j){
			dCnm_vec(counter) = dCnm(i,j);
			++counter;
		}
		
	}

	counter = 0;
	for (int i = 1; i <= degree; ++i){
		for (int j = 1; j <= i; ++j){
			dSnm_vec(counter) = dSnm(i,j);
			++counter;
		}
		
	}


	spherical_harmonics -> SetCnm(Cnm + dCnm);
	spherical_harmonics -> SetSnm(Snm + dSnm);

	sharm_acc_pert = spherical_harmonics -> GetAcceleration(pos);
	dacc_true = sharm_acc_pert - sharm_acc;

	arma::mat partial_C,partial_S;
	spherical_harmonics -> GetPartialHarmonics(pos,partial_C, partial_S);

	dacc_lin = partial_C * dCnm_vec + partial_S * dSnm_vec;


	assert(arma::norm(dacc_lin - dacc_true ) / arma::norm(dacc_true) * 100 < 1e-2);











	
	std::cout << "- Done running test_spherical_harmonics_partials_consistency" << std::endl;

}


void TestsSBCore::test_frame_conversion(){


	std::cout << "- Running test_frame_conversion ..." << std::endl;


	SBGATFrameGraph frame_graph;

	frame_graph.add_frame("N");
	frame_graph.add_frame("B");
	frame_graph.add_transform("N","B");
	arma::vec::fixed<3> mrp = {0.1,0.2,0.3};
	arma::vec::fixed<3> origin = {-0.1,0.2,-0.3};

	frame_graph.set_transform_mrp("N","B",mrp);
	frame_graph.set_transform_origin("N","B",origin);


	arma::vec::fixed<3> zero = {0,0,0};

	assert(arma::norm(frame_graph.convert(zero,"B","N") - origin) == 0);
	
	// TODO : complete with rotation tests
	std::cout << "- Done running test_frame_conversion" << std::endl;

}

void TestsSBCore::test_sbgat_shape_uq(){
	std::cout << "- Running test_sbgat_shape_uq ..." << std::endl;
	
	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName("../../resources/shape_models/skewed.obj");
	reader -> Update(); 
	vtkSmartPointer<vtkPolyData> shape = reader -> GetOutput();

	SBGATTransformShape::ShiftShapeToBarycenter(shape);

	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	mass_properties -> SetInputData(shape);
	mass_properties -> Update();

	std::cout << "Volume: \n" << mass_properties -> GetVolume();
	std::cout << "\nCenter-of-mass: \n" << mass_properties -> GetCenterOfMass().t();
	std::cout << "\nUnit-density inertia: \n" << mass_properties -> GetUnitDensityInertiaTensor();
	std::cout << "\nPrincipal axes MRP: \n" << RBK::dcm_to_mrp(mass_properties -> GetPrincipalAxes()).t();

	vtkSmartPointer<SBGATShapeUncertainty> shape_uq = vtkSmartPointer<SBGATShapeUncertainty>::New();
	shape_uq -> SetInputData(shape);
	shape_uq -> Update();

	auto start = std::chrono::system_clock::now();
	// 
	shape_uq -> ComputeInertiaStatistics(1,0.01);
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "-- Done computing analytical shape uncertainty in " << elapsed_seconds.count() << " s\n";

	double volume_variance_lin = shape_uq -> GetVolumeVariance();
	arma::mat::fixed<3,3> com_covariance_lin = shape_uq -> GetCOMCovariance();
	arma::mat::fixed<6,6> inertia_covariance_lin = shape_uq -> GetInertiaParametrizationCovariance();
	arma::mat::fixed<3,3> principal_dims_covariance_lin = shape_uq -> GetPrincipalDimensionsCovariance();
	arma::mat::fixed<4,4> principal_moments_covariance_lin = shape_uq -> GetPrincipalMomentsCovariance();
	arma::mat::fixed<3,3> mrp_covariance_lin = shape_uq -> GetPrincipalAxesMRPCovariance();

	std::cout << "\nVolume variance: \n";
	std::cout << volume_variance_lin << std::endl;

	std::cout << "\nCOM covariance: \n";
	std::cout << com_covariance_lin << std::endl;

	std::cout << "\nUnit-density Inertia parametrization covariance: \n";
	std::cout << inertia_covariance_lin << std::endl;

	std::cout << "\nPrincipal dimensions covariance: \n";
	std::cout << principal_dims_covariance_lin << std::endl;

	std::cout << "\nPrincipal moments covariance: \n";
	std::cout << principal_moments_covariance_lin << std::endl;

	std::cout << "\nMRP covariance: \n";
	std::cout << mrp_covariance_lin << std::endl;

	start = std::chrono::system_clock::now();
	shape_uq -> ComputeInertiaStatisticsMC(100,1,0.01);
	end = std::chrono::system_clock::now();
	
	elapsed_seconds = end-start;

	std::cout << "\n-- Error from linearized volume variance after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << (shape_uq -> GetVolumeVariance() - volume_variance_lin)/volume_variance_lin * 100 << " %\n";
	std::cout << "-- Max Error from linearized center-of-mass covariance after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetCOMCovariance() - com_covariance_lin).max()/arma::abs(com_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized inertia parametrization covariance after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetInertiaParametrizationCovariance() - inertia_covariance_lin).max()/arma::abs(inertia_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal moments covariance after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalMomentsCovariance() - principal_moments_covariance_lin).max()/arma::abs(principal_moments_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal dimensions covariance after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalDimensionsCovariance() - principal_dims_covariance_lin).max()/arma::abs(principal_dims_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized mrp covariances after 100 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalAxesMRPCovariance() - mrp_covariance_lin).max()/arma::abs(mrp_covariance_lin).max() * 100 << " %\n";

	start = std::chrono::system_clock::now();
	shape_uq -> ComputeInertiaStatisticsMC(1000,1,0.01);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;


	std::cout << "\n-- Error from linearized volume variance after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : "<< (shape_uq -> GetVolumeVariance() - volume_variance_lin)/volume_variance_lin * 100 << " %\n";
	std::cout << "-- Max Error from linearized center-of-mass covariances after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetCOMCovariance() - com_covariance_lin).max()/arma::abs(com_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized inertia parametrization covariance after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetInertiaParametrizationCovariance() - inertia_covariance_lin).max()/arma::abs(inertia_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal moments covariances after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalMomentsCovariance() - principal_moments_covariance_lin).max()/arma::abs(principal_moments_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal dimensions covariance after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalDimensionsCovariance() - principal_dims_covariance_lin).max()/arma::abs(principal_dims_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized mrp covariances after 1000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalAxesMRPCovariance() - mrp_covariance_lin).max()/arma::abs(mrp_covariance_lin).max() * 100 << " %\n";

	start = std::chrono::system_clock::now();
	shape_uq -> ComputeInertiaStatisticsMC(10000,1,0.01);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;

	std::cout << "\n-- Error from linearized volume variance after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : "<< (shape_uq -> GetVolumeVariance() - volume_variance_lin)/volume_variance_lin * 100 << " %\n";
	std::cout << "-- Max Error from linearized center-of-mass covariances after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetCOMCovariance() - com_covariance_lin).max()/arma::abs(com_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized inertia parametrization covariance after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetInertiaParametrizationCovariance() - inertia_covariance_lin).max()/arma::abs(inertia_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal moments covariances after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalMomentsCovariance() - principal_moments_covariance_lin).max()/arma::abs(principal_moments_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal dimensions covariance after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalDimensionsCovariance() - principal_dims_covariance_lin).max()/arma::abs(principal_dims_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized mrp covariances after 10000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalAxesMRPCovariance() - mrp_covariance_lin).max()/arma::abs(mrp_covariance_lin).max() * 100 << " %\n";

	start = std::chrono::system_clock::now();
	shape_uq -> ComputeInertiaStatisticsMC(100000,1,0.01);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;

	std::cout << "\n-- Error from linearized volume variance after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : "<< (shape_uq -> GetVolumeVariance() - volume_variance_lin)/volume_variance_lin * 100 << " %\n";
	std::cout << "-- Max Error from linearized center-of-mass covariances after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetCOMCovariance() - com_covariance_lin).max()/arma::abs(com_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized inertia parametrization covariance after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetInertiaParametrizationCovariance() - inertia_covariance_lin).max()/arma::abs(inertia_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal moments covariances after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalMomentsCovariance() - principal_moments_covariance_lin).max()/arma::abs(principal_moments_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized principal dimensions covariance after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalDimensionsCovariance() - principal_dims_covariance_lin).max()/arma::abs(principal_dims_covariance_lin).max() * 100 << " %\n";
	std::cout << "-- Max Error from linearized mrp covariances after 100000 MC outcomes (in " << elapsed_seconds.count() <<  " seconds) : " << arma::abs(shape_uq -> GetPrincipalAxesMRPCovariance() - mrp_covariance_lin).max()/arma::abs(mrp_covariance_lin).max() * 100 << " %\n";

	std::cout << "\n- Done running test_sbgat_shape_uq ..." << std::endl;

}

void TestsSBCore::test_sbgat_transform_shape(){

	std::cout << "- Running test_sbgat_transform_shape ..." << std::endl;

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

	vtkSmartPointer<vtkPolyData> cube = cleanPolyData -> GetOutput();

	arma::vec x = {1,0,0};
	SBGATTransformShape::Translate(x,cube);

	vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();
	mass_filter -> SetInputData(cube);
	mass_filter -> Update();
	auto com_sbgat = mass_filter -> GetCenterOfMass();

	assert(arma::norm(x - com_sbgat)/arma::norm(x) * 100 < 1e-6);

	SBGATTransformShape::ShiftShapeToBarycenter(cube);
	mass_filter -> SetInputData(cube);
	mass_filter -> Update();

	com_sbgat = mass_filter -> GetCenterOfMass();

	std::cout << "- Done running test_sbgat_transform_shape ..." << std::endl;

}

void TestsSBCore::test_PGM_UQ_partials(){
	std::cout << "- Running test_PGM_UQ_partials ..." << std::endl;
	std::cout << "\t-- Testing ../../resources/shape_models/cube.obj ..." << std::endl;
	std::string filename  = "../../resources/shape_models/cube.obj";
	SBGATPolyhedronGravityModelUQ::TestPartials(filename,5e-2);
	filename  = "../../resources/shape_models/skewed.obj";
	std::cout << "\t-- Testing ../../resources/shape_models/skewed.obj ..." << std::endl;
	SBGATPolyhedronGravityModelUQ::TestPartials(filename,5e-2);
	std::cout << "- Done running test_PGM_UQ_partials ..." << std::endl;
}

void TestsSBCore::test_PGM_UQ_cube(){

	std::cout << "- Running test_PGM_UQ_cube ..." << std::endl;

	int N = 10000;

	arma::arma_rng::set_seed(N) ;

	arma::vec pos = {200,300,400};
	double density = 1970;

	std::string filename  = "../../resources/shape_models/skewed.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(density); 
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	arma::vec::fixed<3> nom_acc;
	double nom_pot;
	pgm_filter -> GetPotentialAcceleration(pos,nom_pot,nom_acc);
	std::cout << "Nominal potential : " << nom_pot << std::endl;
	std::cout << "Nominal acceleration : " << nom_acc.t();


	int N_C = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints();
	arma::mat P_CC = std::pow(1e-1,2) * arma::diagmat<arma::mat>( arma::ones<arma::vec>(3 * N_C));

	arma::vec U_mc(N);
	arma::mat A_mc(3,N);
	arma::mat deviations(3 * N_C,N);

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	for (int i = 0; i < N_C; ++i){
		for (int j = 0; j <= i; ++j){
			const arma::mat::fixed<3,3> & P = P_CC.submat(3 * i,3 * j, 3 * i + 2,3 * j + 2);
			shape_uq.SetCovarianceComponent(P,i,j);
			shape_uq.SetCovarianceComponent(P.t(),j,i);
		}
	}

	arma::mat C_CC = shape_uq.GetCovarianceSquareRoot();
	
	auto start = std::chrono::system_clock::now();	
	double variance_U_analytical = shape_uq.GetVariancePotential(pos);
	arma::mat::fixed<3,3> covariance_A_analytical = shape_uq.GetCovarianceAcceleration(pos);
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Analytical statistics in potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";
	std::cout << "\tAnalytical variance in potential: " << variance_U_analytical << std::endl;
	std::cout << "\tAnalytical covariance in acceleration: \n" << covariance_A_analytical << std::endl;

	// MC
	start = std::chrono::system_clock::now();	
	boost::progress_display progress(N);
	#pragma omp parallel for

	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader_mc = vtkSmartPointer<vtkOBJReader>::New();
		reader_mc -> SetFileName(filename.c_str());
		reader_mc -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner_mc = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner_mc -> SetInputConnection (reader_mc -> GetOutputPort());
		cleaner_mc -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

		cleaner_mc -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_mc = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter_mc -> SetInputConnection(cleaner_mc -> GetOutputPort());
		pgm_filter_mc -> SetDensity(density); 
		pgm_filter_mc -> SetScaleMeters();
		pgm_filter_mc -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq_mc;
		shape_uq_mc.SetPGM(pgm_filter_mc);
		
		deviations.col(i) = C_CC * arma::randn<arma::vec>(3 * N_C);
		shape_uq_mc.ApplyDeviation(deviations.col(i));

		arma::vec::fixed<3> acc;
		double pot;
		shape_uq_mc.GetPGM() -> GetPotentialAcceleration(pos,pot,acc);


		U_mc(i) = pot;
		A_mc.col(i) = acc;

		if (i < 20){
			vtkSmartPointer<SBGATObjWriter> writer = SBGATObjWriter::New();
			writer -> SetInputData(vtkPolyData::SafeDownCast(pgm_filter_mc -> GetInput()));
			std::string path = "../output/cube_mc_shape_" + std::to_string(i) + ".obj";
			writer -> SetFileName(path.c_str());
			writer -> Update();
		}
		++progress;
	}
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	std::cout << "\nMC of potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";

	std::vector<int> steps = {10,100,1000,N};

	for (auto step : steps){

		std::cout << "\t After " << step << " MC outcomes:\n";

		std::cout << "\t\tMC variance in potential: " << arma::var(U_mc.subvec(0,step - 1)) << std::endl;
		std::cout << "\t\tError (%): " << (arma::var(U_mc.subvec(0,step - 1)) - variance_U_analytical)/variance_U_analytical * 100 << std::endl;

		std::cout << "\t\tMC Covariance in acceleration: \n" << arma::cov(A_mc.cols(0,step - 1).t()) << std::endl;
		std::cout << "\t\tError (%): \n" << (arma::cov(A_mc.cols(0,step - 1).t()) - covariance_A_analytical)/covariance_A_analytical * 100 << std::endl;

	}

	
	
	arma::mat P_CC_mc = arma::cov(deviations.t(),deviations.t());
	arma::vec sigma_sq_mc = arma::eig_sym(P_CC_mc);
	arma::vec sigma_sq_true = arma::eig_sym(P_CC);
	
	std::cout << "\t Relative error in MC vs prescribed deviation covariance eigenvalues (%): " << arma::norm(sigma_sq_mc - sigma_sq_true)/arma::norm(sigma_sq_true) * 100 << std::endl;
	std::cout << "\t Absolute difference in MC vs prescribed deviation covariance divided by mean prescribed deviation covariance eigenvalues (%): " << arma::abs(P_CC - P_CC_mc).max()/arma::mean(sigma_sq_true) * 100 << std::endl;


	std::cout << "- Done running test_PGM_UQ_cube ..." << std::endl;


}



void TestsSBCore::test_PGM_UQ_itokawa_m(){

	std::cout << "- Running test_PGM_UQ_itokawa_m ..." << std::endl;

	int N = 10000;
	arma::arma_rng::set_seed(N) ;

	arma::vec pos = {200,300,400};
	double density = 1970;


	std::string filename  = "../../resources/shape_models/itokawa_8_scaled.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(density); 
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	arma::vec::fixed<3> nom_acc;
	double nom_pot;
	pgm_filter -> GetPotentialAcceleration(pos,nom_pot,nom_acc);
	std::cout << "Nominal potential : " << nom_pot << std::endl;
	std::cout << "Nominal acceleration : " << nom_acc.t();


	int N_C = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints();

	arma::vec U_mc(N);
	arma::mat A_mc(3,N);
	arma::mat deviations(3 * N_C,N);

	arma::mat P_CC = std::pow(10e0,2) * arma::diagmat<arma::mat>( arma::ones<arma::vec>(3 * N_C));

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	for (int i = 0; i < N_C; ++i){
		for (int j = 0; j <= i; ++j){
			const arma::mat::fixed<3,3> & P = P_CC.submat(3 * i,3 * j, 3 * i + 2,3 * j + 2);
			shape_uq.SetCovarianceComponent(P,i,j);
			shape_uq.SetCovarianceComponent(P.t(),j,i);
		}
	}

	arma::mat C_CC = shape_uq.GetCovarianceSquareRoot();

	assert(arma::abs(P_CC - C_CC * C_CC.t()).max() < 1e-10);


	auto start = std::chrono::system_clock::now();	
	double variance_U_analytical = shape_uq.GetVariancePotential(pos);
	arma::mat::fixed<3,3> covariance_A_analytical = shape_uq.GetCovarianceAcceleration(pos);
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Analytical statistics in potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";
	std::cout << "\tAnalytical variance in potential: " << variance_U_analytical << std::endl;
	std::cout << "\tAnalytical covariance in acceleration: \n" << covariance_A_analytical << std::endl;

	// MC
	start = std::chrono::system_clock::now();	
	boost::progress_display progress(N);

	#pragma omp parallel for
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader_mc = vtkSmartPointer<vtkOBJReader>::New();
		reader_mc -> SetFileName(filename.c_str());
		reader_mc -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner_mc = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner_mc -> SetInputConnection (reader_mc -> GetOutputPort());
		cleaner_mc -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

		cleaner_mc -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_mc = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter_mc -> SetInputConnection(cleaner_mc -> GetOutputPort());
		pgm_filter_mc -> SetDensity(density); 
		pgm_filter_mc -> SetScaleMeters();
		pgm_filter_mc -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq_mc;
		shape_uq_mc.SetPGM(pgm_filter_mc);
		
		deviations.col(i) = C_CC * arma::randn<arma::vec>(3 * N_C);
		shape_uq_mc.ApplyDeviation(deviations.col(i));

		arma::vec::fixed<3> acc;
		double pot;
		shape_uq_mc.GetPGM() -> GetPotentialAcceleration(pos,pot,acc);

		U_mc(i) = pot;
		A_mc.col(i) = acc;

		if (i < 20){
			vtkSmartPointer<SBGATObjWriter> writer = SBGATObjWriter::New();
			writer -> SetInputData(vtkPolyData::SafeDownCast(pgm_filter_mc -> GetInput()));
			std::string path = "../output/itokawa_m_mc_shape_" + std::to_string(i) + ".obj";
			writer -> SetFileName(path.c_str());
			writer -> Update();
		}
		++progress;
	}
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	std::cout << "\nMC of potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";


	std::vector<int> steps = {10,100,1000,N};

	for (auto step : steps){

		std::cout << "\t After " << step << " MC outcomes:\n";

		std::cout << "\t\tMC variance in potential: " << arma::var(U_mc.subvec(0,step - 1)) << std::endl;
		std::cout << "\t\tError (%): " << (arma::var(U_mc.subvec(0,step - 1)) - variance_U_analytical)/variance_U_analytical * 100 << std::endl;

		std::cout << "\t\tMC Covariance in acceleration: \n" << arma::cov(A_mc.cols(0,step - 1).t()) << std::endl;
		std::cout << "\t\tError (%): \n" << (arma::cov(A_mc.cols(0,step - 1).t()) - covariance_A_analytical)/covariance_A_analytical * 100 << std::endl;

	}
	
	arma::mat P_CC_mc = arma::cov(deviations.t(),deviations.t());
	arma::vec sigma_sq_mc = arma::eig_sym(P_CC_mc);;
	arma::vec sigma_sq_true = arma::eig_sym(P_CC);
	
	std::cout << "\t Relative error in MC vs prescribed deviation covariance eigenvalues (%): " << arma::norm(sigma_sq_mc - sigma_sq_true)/arma::norm(sigma_sq_true) * 100 << std::endl;
	std::cout << "\t Absolute difference in MC vs prescribed deviation covariance divided by mean prescribed deviation covariance eigenvalues (%): " << arma::abs(P_CC - P_CC_mc).max()/arma::mean(sigma_sq_true) * 100 << std::endl;


	std::cout << "- Done running test_PGM_UQ_itokawa_m ..." << std::endl;



}




void TestsSBCore::test_PGM_UQ_itokawa_km(){

	std::cout << "- Running test_PGM_UQ_itokawa_km ..." << std::endl;

	int N = 10000;
	arma::arma_rng::set_seed(N) ;


	arma::vec pos = {200,300,400};
	double density = 1970;


	std::string filename  = "../../resources/shape_models/itokawa_8.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(density); 
	pgm_filter -> SetScaleKiloMeters();
	pgm_filter -> Update();

	arma::vec::fixed<3> nom_acc;
	double nom_pot;

	pgm_filter -> GetPotentialAcceleration(pos,nom_pot,nom_acc);
	std::cout << "Nominal potential : " << nom_pot << std::endl;
	std::cout << "Nominal acceleration : " << nom_acc.t();

	int N_C = vtkPolyData::SafeDownCast(pgm_filter -> GetInput()) -> GetNumberOfPoints();

	arma::vec U_mc(N);
	arma::mat A_mc(3,N);
	arma::mat deviations(3 * N_C,N);

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	shape_uq.ComputeVerticesCovarianceGlobal(10,70);
	arma::mat C_CC_cholesky = shape_uq.GetCovarianceSquareRoot();
	arma::mat C_CC_spectral = shape_uq.GetCovarianceSquareRoot(false);
	arma::mat C_CC;
	arma::mat P_CC = shape_uq.GetVerticesCovariance();

	double error_cholesky = arma::abs(P_CC - C_CC_cholesky * C_CC_cholesky.t()).max();
	double error_spectral = arma::abs(P_CC - C_CC_spectral * C_CC_spectral.t()).max();

	std::cout << "Absolute Error of covariance matrix square root extraction: \n";
	std:cout << "\tCholesky: " << error_cholesky << std::endl;
	std:cout << "\tSpectral decomposition: " << error_spectral << std::endl;

	if (error_cholesky < error_spectral){
		std::cout << "Using cholesky square root";
		C_CC = C_CC_cholesky;
	}
	else{
		std::cout << "Using spectral decomposition square root";
		C_CC = C_CC_spectral;
	}


	auto start = std::chrono::system_clock::now();	
	double variance_U_analytical = shape_uq.GetVariancePotential(pos);
	arma::mat::fixed<3,3> covariance_A_analytical = shape_uq.GetCovarianceAcceleration(pos);
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Analytical statistics in potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";
	std::cout << "\tAnalytical variance in potential: " << variance_U_analytical << std::endl;
	std::cout << "\tAnalytical covariance in acceleration: \n" << covariance_A_analytical << std::endl;

	// MC
	start = std::chrono::system_clock::now();	
	boost::progress_display progress(N);

	#pragma omp parallel for
	for (int i = 0; i < N ; ++i){

		// Reading
		vtkSmartPointer<vtkOBJReader> reader_mc = vtkSmartPointer<vtkOBJReader>::New();
		reader_mc -> SetFileName(filename.c_str());
		reader_mc -> Update(); 

	// Cleaning
		vtkSmartPointer<vtkCleanPolyData> cleaner_mc = vtkSmartPointer<vtkCleanPolyData>::New();
		cleaner_mc -> SetInputConnection (reader_mc -> GetOutputPort());
		cleaner_mc -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

		cleaner_mc -> Update();

	// Creating the PGM dyads
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_mc = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
		pgm_filter_mc -> SetInputConnection(cleaner_mc -> GetOutputPort());
		pgm_filter_mc -> SetDensity(density); 
		pgm_filter_mc -> SetScaleKiloMeters();
		pgm_filter_mc -> Update();

		SBGATPolyhedronGravityModelUQ shape_uq_mc;
		shape_uq_mc.SetPGM(pgm_filter_mc);
		
		deviations.col(i) = C_CC * arma::randn<arma::vec>(3 * N_C);
		shape_uq_mc.ApplyDeviation(deviations.col(i));

		arma::vec::fixed<3> acc;
		double pot;
		shape_uq_mc.GetPGM() -> GetPotentialAcceleration(pos,pot,acc);

		U_mc(i) = pot;
		A_mc.col(i) = acc;

		if (i < 20){
			vtkSmartPointer<SBGATObjWriter> writer = SBGATObjWriter::New();
			writer -> SetInputData(vtkPolyData::SafeDownCast(pgm_filter_mc -> GetInput()));
			std::string path = "../output/itokawa_km_mc_shape_" + std::to_string(i) + ".obj";
			writer -> SetFileName(path.c_str());
			writer -> Update();
		}
		++progress;
	}
	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	std::cout << "\nMC of potential and acceleration computed in " << elapsed_seconds.count() << " seconds\n";


	std::vector<int> steps = {10,100,1000,N};

	for (auto step : steps){

		std::cout << "\t After " << step << " MC outcomes:\n";

		std::cout << "\t\tMC variance in potential: " << arma::var(U_mc.subvec(0,step - 1)) << std::endl;
		std::cout << "\t\tError (%): " << (arma::var(U_mc.subvec(0,step - 1)) - variance_U_analytical)/variance_U_analytical * 100 << std::endl;

		std::cout << "\t\tMC Covariance in acceleration: \n" << arma::cov(A_mc.cols(0,step - 1).t()) << std::endl;
		std::cout << "\t\tError (%): \n" << (arma::cov(A_mc.cols(0,step - 1).t()) - covariance_A_analytical)/covariance_A_analytical * 100 << std::endl;

	}
	
	arma::mat P_CC_mc = arma::cov(deviations.t(),deviations.t());
	arma::vec sigma_sq_mc = arma::eig_sym(P_CC_mc);
	arma::vec sigma_sq_true = arma::eig_sym(P_CC);
	
	
	std::cout << "\t Relative error in MC vs prescribed deviation covariance eigenvalues (%): " << arma::norm(sigma_sq_mc - sigma_sq_true)/arma::norm(sigma_sq_true) * 100 << std::endl;
	std::cout << "\t Absolute difference in MC vs prescribed deviation covariance divided by mean prescribed deviation covariance eigenvalues (%): " << arma::abs(P_CC - P_CC_mc).max()/arma::mean(sigma_sq_true) * 100 << std::endl;


	std::cout << "- Done running test_PGM_UQ_itokawa_km ..." << std::endl;



}


void TestsSBCore::test_PGM_UQ_covariance_consistency(){

	std::cout << "- Running test_PGM_UQ_itokawa_km ..." << std::endl;




	double density = 1970;

	std::string filename  = "../../resources/shape_models/itokawa_8.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Cleaning
	vtkSmartPointer<vtkCleanPolyData> cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
	cleaner -> SetInputConnection (reader -> GetOutputPort());
	cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

	cleaner -> Update();

	// Creating the PGM dyads
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(cleaner -> GetOutputPort());
	pgm_filter -> SetDensity(density); 
	pgm_filter -> SetScaleKiloMeters();
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ shape_uq;
	shape_uq.SetPGM(pgm_filter);

	shape_uq.ComputeVerticesCovarianceGlobal(10,50);

	shape_uq.SaveNonZeroVerticesCovariance("../output/shape_covariance.json");
	shape_uq.SaveVerticesCovariance("../output/shape_covariance.txt");

	arma::mat P_CC;
	P_CC = shape_uq.GetVerticesCovariance();

	if (P_CC.n_rows == 1){
		P_CC = shape_uq.GetVerticesCovariance(false);
	}

	assert(shape_uq.LoadVerticesCovarianceFromJson("../output/shape_covariance.json"));

	arma::mat P_CC_from_Json = shape_uq.GetVerticesCovariance();


	assert(arma::abs(P_CC - P_CC_from_Json).max() / std::pow(10,2) < 1e-10);
	
	std::cout << "- Done running test_PGM_UQ_itokawa_km ..." << std::endl;


}





// TODO : reimplement these tests

// /**
// This test computes simulated radar images for benchmarking purposes
// */
// void TestsSBCore::test_radar_obs(){

// 	std::cout << "- Running test_radar_obs ..." << std::endl;

// 	// Loading in the shape model
// 	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
// 	reader -> SetFileName("../../resources/shape_models/KW4Alpha.obj");
// 	reader -> Update(); 

// 	// Creating the radar object
// 	vtkSmartPointer<SBGATObsRadar> radar = vtkSmartPointer<SBGATObsRadar>::New();
// 	radar -> SetInputConnection(reader -> GetOutputPort());
// 	radar -> SetScaleKiloMeters();
// 	radar -> Update();


// 	// Arguments
// 	arma::vec spin = {0,0,1};
// 	arma::vec dir = {1,0,0};
// 	double period = 4 * 3600; // 4 hours
// 	int images = 48; 
// 	int N = 100;

// 	double r_bin = 7.5;//(m)
// 	double rr_bin = 7.9e-3;//(m/s)

// 	// A sequence of images is collected
// 	SBGATRadarObsSequence measurement_sequence;
// 	auto start = std::chrono::system_clock::now();

// 	for (int i  = 0; i < images; ++i){
// 		double dt = 1.5 * double(i) / ((images - 1)) * period;

// 		std::cout << " --- Ray tracing " +std::to_string(i + 1) + "/" +std::to_string(images) + " ...\n";

// 		radar -> CollectMeasurementsSimpleSpin(measurement_sequence,N,dt,period,dir,spin);
// 	}

// 	radar -> BinObservations(measurement_sequence,r_bin,rr_bin);
// 	std::cout << " --- Done binning ...\n";

// 	radar -> SaveImages("../radar_output/");
// 	std::cout << " --- Done saving ...\n";


// 	auto end = std::chrono::system_clock::now();

// 	std::chrono::duration<double> elapsed_seconds = end-start;
// 	std::cout << "-- Done collecting radar images in " << elapsed_seconds.count() << " s\n";
// 	std::cout << "-- test_radar_obs successful" << std::endl;

// }


// *
// This test computes simulated lightcurves for benchmarking purposes

// void TestsSBCore::test_lightcurve_obs(){

// 	std::cout << "- Running test_lightcurve_obs ..." << std::endl;

// 	// Loading in the shape model
// 	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
// 	reader -> SetFileName("../../resources/shape_models/KW4Alpha.obj");
// 	reader -> Update(); 

// 	// Creating the radar object
// 	vtkSmartPointer<SBGATObsLightcurve> lightcurve = vtkSmartPointer<SBGATObsLightcurve>::New();
// 	lightcurve -> SetInputConnection(reader -> GetOutputPort());
// 	lightcurve -> SetScaleKiloMeters();
// 	lightcurve -> Update();


// 	// Arguments
// 	arma::vec spin = {0,0,1};
// 	arma::vec target_pos = {1e6,0,0};
// 	arma::vec observer_pos = {1e6,1e6,0};
// 	double period = 4 * 3600; // 4 hours
// 	int images = 100; 
// 	int N = 100;

// 	// A sequence of images is collected
// 	std::vector<std::array<double, 2> > measurements;
// 	auto start = std::chrono::system_clock::now();

// 	for (int i  = 0; i < images; ++i){
// 		double dt = 360 * double(i) ; // one data point evey 360 s

// 		std::cout << " --- Ray tracing " +std::to_string(i + 1) + "/" +std::to_string(images) + " ...\n";
// 		lightcurve -> CollectMeasurementsSimpleSpin(measurements,N,dt,period,target_pos,observer_pos,spin);
// 		std::cout << measurements.back()[1] << std::endl;
// 	}


// 	auto end = std::chrono::system_clock::now();

// 	std::chrono::duration<double> elapsed_seconds = end-start;
// 	std::cout << "-- Done collecting lightcurve " << elapsed_seconds.count() << " s\n";
// 	std::cout << "-- test_lightcurve_obs successful" << std::endl;

// 	arma::vec mes_arma(measurements.size());
// 	for (unsigned int i = 0; i < measurements.size(); ++i){
// 		mes_arma(i) = measurements[i][1];
// 	}
// 	mes_arma.save("lightcurve.txt",arma::raw_ascii);

// }












