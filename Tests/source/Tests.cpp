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
	TestsSBCore::test_spherical_harmonics_partials_consistency();
	TestsSBCore::test_frame_conversion();
	TestsSBCore::test_sbgat_mass_properties();
	TestsSBCore::test_sbgat_pgm();
	TestsSBCore::test_sbgat_pgm_speed();
	TestsSBCore::test_spherical_harmonics_coefs_consistency();

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
	auto inertia_moments_sbgat = mass_filter -> GetInertiaMoments();

	assert(arma::norm(inertia_moments - inertia_moments_sbgat)/arma::norm(inertia_moments_sbgat) < 1e-8);
	assert(arma::norm(translation_vector - com_sbgat)/arma::norm(com_sbgat) < 1e-8);

	std::cout << "- Done running test_sbgat_mass_properties" << std::endl;

}


/**
This test computes the surface accelerations at the center of each facet over a polydata
of Eros for benchmarking purposes
*/
void TestsSBCore::test_sbgat_pgm_speed(){
	std::cout << "- Running test_sbgat_pgm_speed ..." << std::endl;
	
	std::string filename  = "../input/KW4Alpha.obj";

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(filename.c_str());
	reader -> Update(); 

	// Creating the PGM dyads
	std::cout << "-- Creating dyads...\n";
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(1970); // density in kg/m^3
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
	arma::mat surface_accelerations(cellCentersFilter -> GetOutput() -> GetNumberOfPoints(),3);

	assert(polydata -> GetNumberOfCells() == cellCentersFilter -> GetOutput() -> GetNumberOfPoints());


	auto start = std::chrono::system_clock::now();
	std::cout << "-- Computing pgm accelerations at " << cellCentersFilter -> GetOutput() -> GetNumberOfPoints() << " facet centers over the surface of " << filename <<  " . This may take a few minutes ...\n";
	
	boost::progress_display progress(cellCentersFilter -> GetOutput() -> GetNumberOfPoints());

	for (vtkIdType i = 0; i < cellCentersFilter -> GetOutput() -> GetNumberOfPoints(); ++i){
		double p[3];
		cellCentersFilter -> GetOutput() -> GetPoint(i, p);
		surface_accelerations.row(i) = pgm_filter -> GetAcceleration(p ).t();
		++progress;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "-- Done computing pgm accelerations in " << elapsed_seconds.count() << " s\n";


	// The surface acceleration are re-calculated, this time using the static method SBGATPolyhedronGravityModel::ComputeSurfacePGM
	// that will also evaluate the inertial potential and geopotential in the process
	std::vector<unsigned int> queried_elements;
	for (int i = 0; i < polydata -> GetNumberOfCells(); ++i){
		queried_elements.push_back(i);
	}

	// Rotation period of KW4
	arma::vec omega = {0,0,0.0006312};

	std::vector<double> slopes,inertial_potentials,
	body_fixed_potentials,
	inertial_acc_magnitudes,
	body_fixed_acc_magnitudes;

	SBGATPolyhedronGravityModel::ComputeSurfacePGM(polydata,
		queried_elements,
		false,
		pgm_filter -> GetDensity(),
		omega,
		slopes,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes);


	SBGATPolyhedronGravityModel::SaveSurfacePGM(polydata,
		queried_elements,
		false,
		pgm_filter -> GetMass(),
		omega,
		slopes,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes,
		"surface_pgm.json");










	std::cout << "- Done running test_sbgat_pgm_speed" << std::endl;

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


		arma::vec::fixed<3> p4_perturbed = (1 + 1e-5) * p4;
		

		arma::vec::fixed<3> pgm_acc2_perturbed = pgm_filter -> GetAcceleration(p4_perturbed);


		auto d_acc_true = pgm_acc2_perturbed - pgm_acc2;
		arma::vec::fixed<3> d_acc_linear = pgm_gravity_gradient_matrix * (p4_perturbed - p4);
		assert(100 * arma::norm(d_acc_linear - d_acc_true)/arma::norm(d_acc_true) < 1e-2);

	}

	std::cout << "- Done running test_sbgat_pgm" << std::endl;

}


/**
This test checks the consistency of the spherical harmonics coefficients computation 
about a shape model of KW4 
*/

void TestsSBCore::test_spherical_harmonics_coefs_consistency() {

	std::cout << "- Running test_spherical_harmonics_coefs_consistency ..." << std::endl;


	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName("../input/KW4Alpha.obj");
	reader -> Update(); 

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of KW4 (kg/m^3). 
	double density = 2000.0;

	// Reference radius of KW4 (km)
	double ref_radius = 1.317/2;


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
	spherical_harmonics -> SaveToJson("../gravity_output/harmo.json");


	// The accelerations are evaluated at the query point
	arma::vec::fixed<3>  pos = {3,5,-2};
	arma::vec::fixed<3> pgm_acc = pgm_filter -> GetAcceleration(pos.colptr(0));
	arma::vec::fixed<3> sharm_acc = spherical_harmonics -> GetAcceleration(pos);


	// The spherical harmonics are read from the just-saved JSON file and re-evaluated
	vtkSmartPointer<SBGATSphericalHarmo> spherical_harmonics_from_file = vtkSmartPointer<SBGATSphericalHarmo>::New();
	spherical_harmonics_from_file -> LoadFromJson("../gravity_output/harmo.json");
	arma::vec::fixed<3> sharm_acc_from_file = spherical_harmonics_from_file -> GetAcceleration(pos);

	// The accelerations should be consistent with the one previously computed
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
	reader -> SetFileName("../input/KW4Alpha.obj");
	reader -> Update(); 

	// Harmonics up to degree five are computed
	int degree = 5;

	// Density of KW4 (kg/m^3). 
	double density = 2000.0;

	// Reference radius of KW4 (km)
	double ref_radius = 1.317/2;


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
	arma::vec::fixed<3> pos = {3,5,-2};
	arma::vec::fixed<3> dpos = 1e-3 * pos;

	arma::vec::fixed<3> sharm_acc = spherical_harmonics -> GetAcceleration(pos);
	arma::vec::fixed<3> sharm_acc_pert = spherical_harmonics -> GetAcceleration(pos + dpos);
	arma::vec::fixed<3> dacc_true = sharm_acc_pert - sharm_acc;
	arma::mat::fixed<3,3> gravity_gradient_mat;
	spherical_harmonics -> GetGravityGradientMatrix(pos,gravity_gradient_mat);

	arma::vec::fixed<3> dacc_lin = gravity_gradient_mat * dpos;
	assert(arma::norm(dacc_lin - dacc_true ) / arma::norm(dacc_true) * 100 < 2e-1);


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






// TODO : reimplement these tests

// /**
// This test computes simulated radar images for benchmarking purposes
// */
// void TestsSBCore::test_radar_obs(){

// 	std::cout << "- Running test_radar_obs ..." << std::endl;

// 	// Loading in the shape model
// 	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
// 	reader -> SetFileName("../input/KW4Alpha.obj");
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
// 	reader -> SetFileName("../input/KW4Alpha.obj");
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












