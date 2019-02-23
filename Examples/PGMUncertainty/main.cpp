#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>

#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>


#include <json.hpp>


int main(){

	std::ifstream i("input_file.json");
	nlohmann::json input_data;
	i >> input_data;

	std::string PATH_SHAPE = input_data["PATH_SHAPE"];
	double CORRELATION_DISTANCE =  input_data["CORRELATION_DISTANCE"];

	double ERROR_STANDARD_DEV  = input_data["ERROR_STANDARD_DEV"];
	double DENSITY  = input_data["DENSITY"];
	bool UNIT_IN_METERS  = input_data["UNIT_IN_METERS"];

	int N_monte_carlo = input_data["N_MONTE_CARLO"];
	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];

	std::cout << "Path to shape: " << PATH_SHAPE << std::endl;
	std::cout << "Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "Monte Carlo Draws : " << N_monte_carlo << std::endl;

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
	pgm_filter -> Update();

	SBGATPolyhedronGravityModelUQ pgm_uq;
	pgm_uq.SetPGM(pgm_filter);

	// Populate the covariances
	pgm_uq.ComputeCovariancesGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);

	// Save the covariances
	pgm_uq.SaveNonZeroCovariances(OUTPUT_DIR + "shape_covariance.json");



	return 0;
}
