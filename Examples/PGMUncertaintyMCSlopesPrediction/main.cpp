#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>
#include <SBGATObjWriter.hpp>

#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>

#include <json.hpp>
#include <boost/progress.hpp>

#include <armadillo>

int main(){


	std::ifstream i("input_file.json");
	nlohmann::json input_data;
	i >> input_data;

	arma::arma_rng::set_seed(0);

	std::string PATH_SHAPE = input_data["PATH_SHAPE"];
	double CORRELATION_DISTANCE =  input_data["CORRELATION_DISTANCE"];

	double ERROR_STANDARD_DEV  = input_data["ERROR_STANDARD_DEV"];
	double DENSITY  = input_data["DENSITY"];
	double PERIOD_SD  = input_data["PERIOD_SD"];
	double PERIOD  = input_data["PERIOD"];
	
	bool UNIT_IN_METERS  = input_data["UNIT_IN_METERS"];
	bool HOLD_MASS_CONSTANT  = input_data["HOLD_MASS_CONSTANT"];


	int N_MONTE_CARLO = input_data["N_MONTE_CARLO"];
	std::vector<int> COV_REGION_CENTERS = input_data["COV_REGION_CENTERS"];

	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];
	std::string UNCERTAINTY_TYPE = input_data["UNCERTAINTY_TYPE"];



	std::cout << "- Path to shape: " << PATH_SHAPE << std::endl;
	std::cout << "- Uncertainty type: " << UNCERTAINTY_TYPE << std::endl;
	std::cout << "- Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "- Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "- Standard deviation on rotation period (s) : " << PERIOD_SD << std::endl;
	std::cout << "- Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "- Rotation period (s) : " << PERIOD << std::endl;
	std::cout << "- Covariance region centers:\n" ;
	for(auto center : COV_REGION_CENTERS){
		std::cout << "\t" << center << std::endl;
	}
	

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(PATH_SHAPE.c_str());
	reader -> Update(); 

	// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(DENSITY);

	
	if(UNIT_IN_METERS){
		pgm_filter -> SetScaleMeters();
	} else{
		pgm_filter -> SetScaleKiloMeters();
	}

	std::cout << "Building pgm ...\n";
	pgm_filter -> SetOmega(2 * arma::datum::pi / PERIOD * arma::vec({0,0,1}));
	pgm_filter -> Update();

	double mass = DENSITY * pgm_filter -> GetVolume();
	arma::vec::fixed<3> omega = pgm_filter -> GetOmega();
	vtkSmartPointer<vtkPolyData> polydata = reader -> GetOutput();

	// An instance of SBGATPolyhedronGravityModelUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	SBGATPolyhedronGravityModelUQ pgm_uq;
	pgm_uq.SetModel(pgm_filter);
	pgm_uq.SetPeriodErrorStandardDeviation(PERIOD_SD);

	pgm_uq.PrecomputeMassPropertiesPartials();

	// Saving baseline slices
	pgm_uq.TakeAndSaveSlice(0,OUTPUT_DIR + "baseline_slice_x.txt",0);
	pgm_uq.TakeAndSaveSlice(1,OUTPUT_DIR + "baseline_slice_y.txt",0);
	pgm_uq.TakeAndSaveSlice(2,OUTPUT_DIR + "baseline_slice_z.txt",0);

	std::cout << "Populating shape covariance ...\n";

	// Populate the shape vertices covariance
	
	if (UNCERTAINTY_TYPE == "radial"){
		for (auto region_center : COV_REGION_CENTERS){
			pgm_uq.AddRadialUncertaintyRegionToCovariance(region_center,ERROR_STANDARD_DEV,CORRELATION_DISTANCE);

		}
	}
	else if (UNCERTAINTY_TYPE == "normal"){
		for (auto region_center : COV_REGION_CENTERS){

			pgm_uq.AddRadialUncertaintyRegionToCovariance(region_center,ERROR_STANDARD_DEV,CORRELATION_DISTANCE);

		}
	}
	else if (UNCERTAINTY_TYPE == "global"){
		pgm_uq.ComputeVerticesCovarianceGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);
	}
	else{
		throw(std::runtime_error("Got unknown uncertainty direction type: " + UNCERTAINTY_TYPE));
	}


	arma::mat P_CC = pgm_uq.GetVerticesCovariance();
	arma::mat C_CC = pgm_uq.GetCovarianceSquareRoot();
	
	P_CC.save(OUTPUT_DIR + "full_covariance.txt",arma::raw_ascii);
	C_CC.save(OUTPUT_DIR + "full_covariance_sqrt.txt",arma::raw_ascii);

	// Regularizing the covariance
	int regularized_eigen_values = pgm_uq.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";
	C_CC = pgm_uq.GetCovarianceSquareRoot();
	
	P_CC = pgm_uq.GetVerticesCovariance();
	P_CC.save(OUTPUT_DIR + "full_covariance_regularized.txt",arma::raw_ascii);
	C_CC.save(OUTPUT_DIR + "full_covariance_sqrt_regularized.txt",arma::raw_ascii);

	regularized_eigen_values = pgm_uq.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";

	std::cout << "Maximum absolute error in covariance square root: " << arma::abs(P_CC - C_CC * C_CC.t()).max() << std::endl;

	std::cout << "Saving non-zero partition of shape covariance ...\n";
	



	std::vector<double> slopes;
	std::vector<double> inertial_potentials;
	std::vector<double> body_fixed_potentials;
	std::vector<double> inertial_acc_magnitudes;
	std::vector<double> body_fixed_acc_magnitudes;


	std::vector<unsigned int> queried_elements;
	for (int i = 0; i < pgm_filter -> GetN_facets() ; ++i){
		queried_elements.push_back(i);
	}


	SBGATPolyhedronGravityModel::ComputeSurfacePGM(polydata,
		queried_elements,
		UNIT_IN_METERS,
		DENSITY,
		omega,
		slopes,
		inertial_potentials,
		body_fixed_potentials,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes);

	// Analytical UQ
	std::vector<double> analytical_variances_slopes;

	std::cout << "Computing analytical uncertainties ... ";
	auto start = std::chrono::system_clock::now();

	pgm_uq.GetVarianceSlopes(analytical_variances_slopes,queried_elements,HOLD_MASS_CONSTANT);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "Done computing analytical uncertainties in " << elapsed_seconds.count() << " s\n";

	std::vector<double> sd_slopes(analytical_variances_slopes.size());
	
	for (unsigned int i = 0; i < analytical_variances_slopes.size(); ++i){
		sd_slopes[i] = std::sqrt(analytical_variances_slopes[i]);
	}



	SBGATPolyhedronGravityModel::SaveSurfacePGM(polydata,
		queried_elements,
		UNIT_IN_METERS,
		mass,
		omega,
		slopes,
		inertial_potentials,
		body_fixed_potentials,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes,
		sd_slopes,
		OUTPUT_DIR + "surface_pgm.json");


	


	return 0;
}
