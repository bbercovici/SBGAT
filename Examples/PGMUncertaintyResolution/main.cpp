#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>
#include <SBGATObjWriter.hpp>

#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>

#include <json.hpp>
#include <boost/progress.hpp>

#include <armadillo>

int main(){

	std::ifstream input("input_file.json");
	nlohmann::json input_data;
	input >> input_data;

	std::string ITOKAWA_8_SHAPE = input_data["ITOKAWA_8_SHAPE"];
	std::string ITOKAWA_32_SHAPE = input_data["ITOKAWA_32_SHAPE"];
	double CORRELATION_DISTANCE =  input_data["CORRELATION_DISTANCE"];
	double ERROR_STANDARD_DEV  = input_data["ERROR_STANDARD_DEV"];
	
	double DENSITY  = input_data["DENSITY"];
	bool UNIT_IN_METERS  = input_data["UNIT_IN_METERS"];
	bool HOLD_MASS_CONSTANT  = input_data["HOLD_MASS_CONSTANT"];
	int PROJECTION_AXIS = input_data["PROJECTION_AXIS"];
	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];
	std::string UNCERTAINTY_TYPE = input_data["UNCERTAINTY_TYPE"];
	
	std::cout << "- Path to ITOKAWA_8_SHAPE: " << ITOKAWA_8_SHAPE << std::endl;
	std::cout << "- Path to ITOKAWA_32_SHAPE: " << ITOKAWA_32_SHAPE << std::endl;

	std::cout << "- Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "- Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "- Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "- Uncertainty type: " << UNCERTAINTY_TYPE << std::endl;
	std::cout << "- Projection axis : " << PROJECTION_AXIS << std::endl;
	


	// Reading
	vtkSmartPointer<vtkOBJReader> reader_8 = vtkSmartPointer<vtkOBJReader>::New();
	reader_8 -> SetFileName(ITOKAWA_8_SHAPE.c_str());
	reader_8 -> Update(); 

	// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_itokawa_8 = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter_itokawa_8 -> SetInputConnection(reader_8 -> GetOutputPort());
	pgm_filter_itokawa_8 -> SetDensity(DENSITY);


	if(UNIT_IN_METERS){
		pgm_filter_itokawa_8 -> SetScaleMeters();

	} else{
		pgm_filter_itokawa_8 -> SetScaleKiloMeters();

	}
	pgm_filter_itokawa_8 -> Update();

	// Reading
	vtkSmartPointer<vtkOBJReader> reader_32 = vtkSmartPointer<vtkOBJReader>::New();

	reader_32 -> SetFileName(ITOKAWA_32_SHAPE.c_str());
	reader_32 -> Update(); 

	// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter_itokawa_32 = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter_itokawa_32 -> SetInputConnection(reader_32 -> GetOutputPort());
	pgm_filter_itokawa_32 -> SetDensity(DENSITY);


	
	if(UNIT_IN_METERS){
		pgm_filter_itokawa_32 -> SetScaleMeters();

	} else{
		pgm_filter_itokawa_32 -> SetScaleKiloMeters();

	}
	pgm_filter_itokawa_32 -> Update();


	// An instance of SBGATPolyhedronGravityModelUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	SBGATPolyhedronGravityModelUQ pgm_uq_itokawa_8;
	pgm_uq_itokawa_8.SetModel(pgm_filter_itokawa_8);
	pgm_uq_itokawa_8.PrecomputeMassPropertiesPartials();
	
	// Populate the shape vertices covariance

	std::cout << "Computing covariance\n";
	pgm_uq_itokawa_8.ComputeVerticesCovarianceGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);
	

	arma::mat P_CC = pgm_uq_itokawa_8.GetVerticesCovariance();
	arma::mat C_CC = pgm_uq_itokawa_8.GetCovarianceSquareRoot();
	std::cout << "Saving covariance\n";
		
	P_CC.save(OUTPUT_DIR + "full_covariance.txt",arma::raw_ascii);
	C_CC.save(OUTPUT_DIR + "full_covariance_sqrt.txt",arma::raw_ascii);
	std::cout << "Regularizing covariance\n";

	// Regularizing the covariance
	int regularized_eigen_values = pgm_uq_itokawa_8.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";
	C_CC = pgm_uq_itokawa_8.GetCovarianceSquareRoot();
	
	P_CC = pgm_uq_itokawa_8.GetVerticesCovariance();
	P_CC.save(OUTPUT_DIR + "full_covariance_regularized.txt",arma::raw_ascii);
	C_CC.save(OUTPUT_DIR + "full_covariance_sqrt_regularized.txt",arma::raw_ascii);

	regularized_eigen_values = pgm_uq_itokawa_8.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";
	std::cout << "Maximum absolute error in covariance square root: " << arma::abs(P_CC - C_CC * C_CC.t()).max() << std::endl;

	// Saving baseline slices
	pgm_uq_itokawa_8.TakeAndSaveSlice(0,OUTPUT_DIR + "baseline_slice_x.txt",1e-6);
	pgm_uq_itokawa_8.TakeAndSaveSlice(1,OUTPUT_DIR + "baseline_slice_y.txt",1e-6);
	pgm_uq_itokawa_8.TakeAndSaveSlice(2,OUTPUT_DIR + "baseline_slice_z.txt",1e-6);

	// Running a Monte Carlo on a subset of the grid

	arma::vec::fixed<3> e0 = {1,0,0};
	arma::vec::fixed<3> e1 = {0,1,0};
	arma::vec::fixed<3> e2 = {0,0,1};
	arma::vec::fixed<3> e3 = arma::normalise(arma::vec({1,1,0}));
	arma::vec::fixed<3> e4 = arma::normalise(arma::vec({0,1,1}));
	arma::vec::fixed<3> e5 = arma::normalise(arma::vec({1,0,1}));
	arma::vec::fixed<3> e6 = arma::normalise(arma::vec({-1,1,0}));
	arma::vec::fixed<3> e7 = arma::normalise(arma::vec({0,-1,1}));
	arma::vec::fixed<3> e8 = arma::normalise(arma::vec({1,0,-1}));

	std::vector<arma::vec::fixed<3> > all_positions;
	std::vector<double> distances = {400,500,600,700};
	for (auto dist : distances){
		all_positions.push_back(dist * e0);
		all_positions.push_back(dist * e1);
		all_positions.push_back(dist * e2);
		all_positions.push_back(dist * e3);
		all_positions.push_back(dist * e4);
		all_positions.push_back(dist * e5);
		all_positions.push_back(dist * e6);
		all_positions.push_back(dist * e7);
		all_positions.push_back(dist * e8);
		all_positions.push_back(- dist * e0);
		all_positions.push_back(- dist * e1);
		all_positions.push_back(- dist * e2);
		all_positions.push_back(- dist * e3);
		all_positions.push_back(- dist * e4);
		all_positions.push_back(- dist * e5);
		all_positions.push_back(- dist * e6);
		all_positions.push_back(- dist * e7);
		all_positions.push_back(- dist * e8);
	}




	// For every location in the grid, evaluate the acceleration arising from itokawa_8 and itokawa_32
	// and analytical uncertainty arising from the prescribed uncertainty in itokawa_8. The difference between
	// the two accelerations should match the predicted uncertainty

	
	// For every point in the grid, evaluate the analytical acceleration covariance and 
	// run a monte carlo to get a statistical covariance to compare against
	arma::mat all_positions_arma = arma::zeros<arma::mat>(3,all_positions.size());
	arma::vec all_positions_arma_normalized_acceleration_difference_sd(all_positions.size());
	boost::progress_display progress(all_positions.size());

	auto start = std::chrono::system_clock::now();

	// #pragma omp parallel for
	for (int p = 0; p < all_positions.size(); ++p){


		const arma::vec::fixed<3> & grid_point = all_positions[p];


		arma::vec::fixed<3> acc_8 = pgm_filter_itokawa_8 -> GetAcceleration(grid_point);
		arma::vec::fixed<3> acc_32 = pgm_filter_itokawa_32 -> GetAcceleration(grid_point);
		arma::vec::fixed<3> acceleration_difference =  acc_32 - acc_8;
		arma::mat::fixed<3,3> itokawa_8_acceleration_covariance = pgm_uq_itokawa_8.GetCovarianceAcceleration(grid_point);


		
		all_positions_arma_normalized_acceleration_difference_sd(p) = std::sqrt(arma::dot(acceleration_difference,arma::inv(itokawa_8_acceleration_covariance) * acceleration_difference));
		all_positions_arma.col(p) = grid_point;

		

		++progress;
	}


	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	
	std::cout << "\n-- Done evaluating over grid in " << elapsed_seconds.count() << " s\n";



	all_positions_arma.save(OUTPUT_DIR + "all_positions_arma.txt",arma::raw_ascii);
	all_positions_arma_normalized_acceleration_difference_sd.save(OUTPUT_DIR + "all_positions_arma_normalized_acceleration_difference_sd.txt",arma::raw_ascii);


	return 0;
}
