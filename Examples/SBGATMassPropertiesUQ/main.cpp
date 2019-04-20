#include <SBGATMassProperties.hpp>
#include <SBGATMassPropertiesUQ.hpp>

#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>

#include <json.hpp>
#include <boost/progress.hpp>
#include <chrono>

int main(){


	// Load the input file
	std::ifstream i("input_file.json");
	nlohmann::json input_data;
	i >> input_data;

	// Path to considered shape
	std::string PATH_SHAPE = input_data["PATH_SHAPE"];

	// Correlation distance (always in m)
	double CORRELATION_DISTANCE =  input_data["CORRELATION_DISTANCE"];

	// Standard deviation in normal error (always in m)
	double ERROR_STANDARD_DEV  = input_data["ERROR_STANDARD_DEV"];

	// Bulk density (always in kg/m^3)
	double DENSITY  = input_data["DENSITY"];

	// True is $PATH_SHAPE points to an .obj shape whose coordinates are expressed in 
	// meters, False otherwise
	bool UNIT_IN_METERS  = input_data["UNIT_IN_METERS"];

	// Number of Monte-Carlo samples
	int N_MONTE_CARLO = input_data["N_MONTE_CARLO"];

	// Path to directory in which the output results must be saved
	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];

	std::cout << "- Path to shape: " << PATH_SHAPE << std::endl;
	std::cout << "- Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "- Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "- Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "- Monte-Carlo draws : " << N_MONTE_CARLO << std::endl;

	// Reading shape model
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(PATH_SHAPE.c_str());
	reader -> Update(); 

	// An instance of SBGATMassProperties is created to evaluate the mass properties
	// of the considered shape
	vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();
	mass_filter -> SetInputConnection(reader -> GetOutputPort());
	mass_filter -> SetDensity(DENSITY);

	// Sets the scaling factor of the input shape
	if(UNIT_IN_METERS){
		mass_filter -> SetScaleMeters();
	} else{
		mass_filter -> SetScaleKiloMeters();
	}

	// Computes the volume, center-of-mass,...
	mass_filter -> Update();



	// An instance of SBGATMassPropertiesUQ is created to perform
	// uncertainty quantification in the inertias of the considered shape
	SBGATMassPropertiesUQ mass_uq;

	// The SBGATMassPropertiesUQ is provided with the reference inertias
	mass_uq.SetModel(mass_filter);

	// Populate the shape vertices covariance
	mass_uq.ComputeVerticesCovarianceGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);
	
	// Regularizing the covariance
	int regularized_eigen_values = mass_uq.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";
	
	arma::mat C_CC = mass_uq.GetCovarianceSquareRoot();
	arma::mat P_CC = mass_uq.GetVerticesCovariance();
	std::cout << "Maximum absolute error in covariance square root: " << arma::abs(P_CC - C_CC * C_CC.t()).max() << std::endl;
	
	std::cout << "Saving shape covariance ...\n";
	P_CC.save(OUTPUT_DIR + "full_covariance.txt",arma::raw_ascii);

	// Saving baseline slices
	mass_uq.TakeAndSaveSlice(0,OUTPUT_DIR + "baseline_slice_x.txt",1e-6);
	mass_uq.TakeAndSaveSlice(1,OUTPUT_DIR + "baseline_slice_y.txt",1e-6);
	mass_uq.TakeAndSaveSlice(2,OUTPUT_DIR + "baseline_slice_z.txt",1e-6);



	// The analytical uncertainties are computed
	auto start = std::chrono::system_clock::now();
	mass_uq.PrecomputeMassPropertiesPartials();
	const arma::mat & partialIPartialC = mass_uq.GetPartialIPartialC();
	const arma::mat & partialComPartialC = mass_uq.GetPartialComPartialC();
	const arma::rowvec & partialVolumePartialC = mass_uq.GetPartialVolumePartialC();
	double sigma_vol = std::sqrt(arma::dot(partialVolumePartialC.t(),P_CC * partialVolumePartialC.t()));
	arma::mat cov_com = partialComPartialC * P_CC * partialComPartialC.t();
	arma::mat cov_I = partialIPartialC * P_CC * partialIPartialC.t();

	std::cout << "Standard deviation on volume (m^3): " << sigma_vol << " \n";
	cov_com.print("Covariance in center-of-mass (m^4): ");
	cov_I.print("Covariance in unit-density inertia tensor (m^10): ");

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;
	std::cout << "Time elapsed computing analytical uncertainties: " << elapsed_seconds.count() << " s" << std::endl;
	

	// The containers storing the MC results are created
	arma::vec all_volumes;
	arma::mat all_com;
	arma::mat all_inertia;
	arma::mat deviations;

	start = std::chrono::system_clock::now();

	// A monte-carlo simulation is run
	SBGATMassPropertiesUQ::RunMCUQVolumeCOMInertia(PATH_SHAPE,
		DENSITY,
		UNIT_IN_METERS,
		C_CC,
		N_MONTE_CARLO,
		OUTPUT_DIR,
		std::min(60,N_MONTE_CARLO),
		deviations,
		all_volumes,
		all_com,
		all_inertia);

	double sigma_vol_MC = arma::stddev(all_volumes);
	arma::mat cov_com_MC = arma::cov(all_com.t(),all_com.t());
	arma::mat cov_I_MC = arma::cov(all_inertia.t(),all_inertia.t());

	std::cout << "Standard deviation on volume (m^3): " << sigma_vol_MC << " \n";
	cov_com_MC.print("Covariance in center-of-mass (m^4): ");
	cov_I_MC.print("Covariance in unit-density inertia tensor (m^10): ");
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	std::cout << "Time elapsed computing Monte Carlo distribution: " << elapsed_seconds.count() << " s" << std::endl;

	std::cout << "Error volume (\%): " << (sigma_vol - sigma_vol_MC)/mass_filter -> GetVolume() << std::endl;
	std::cout << "Error center-of-mass (\%): " << (arma::sqrt(arma::abs(cov_com - cov_com_MC))).max()/arma::norm(mass_filter -> GetCenterOfMass()) << std::endl;
	std::cout << "Error Inertia (\%): " << (arma::sqrt(arma::abs(cov_I - cov_I_MC))).max()/arma::trace(mass_filter -> GetUnitDensityInertiaTensor()) << std::endl;


	std::cout << "##############################################################################\n";










	return 0;
}
