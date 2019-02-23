#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>

#include <vtkCleanPolyData.h>
#include <vtkOBJReader.h>

#include <json.hpp>
#include <boost/progress.hpp>


int main(){

	std::ifstream i("input_file.json");
	nlohmann::json input_data;
	i >> input_data;

	std::string PATH_SHAPE = input_data["PATH_SHAPE"];
	double CORRELATION_DISTANCE =  input_data["CORRELATION_DISTANCE"];

	double ERROR_STANDARD_DEV  = input_data["ERROR_STANDARD_DEV"];
	double DENSITY  = input_data["DENSITY"];
	double STEP_SIZE  = input_data["STEP_SIZE"];

	bool UNIT_IN_METERS  = input_data["UNIT_IN_METERS"];

	int N_monte_carlo = input_data["N_MONTE_CARLO"];
	int PROJECTION_AXIS = input_data["PROJECTION_AXIS"];

	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];

	std::cout << "- Path to shape: " << PATH_SHAPE << std::endl;
	std::cout << "- Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "- Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "- Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "- Monte Carlo Draws : " << N_monte_carlo << std::endl;

	// Reading
	vtkSmartPointer<vtkOBJReader> reader = vtkSmartPointer<vtkOBJReader>::New();
	reader -> SetFileName(PATH_SHAPE.c_str());
	reader -> Update(); 

	// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputConnection(reader -> GetOutputPort());
	pgm_filter -> SetDensity(DENSITY);

	// An instance of SBGATPolyhedronGravityModel is created to get the bounding
	// box of the provided polydata
	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	mass_properties -> SetInputConnection(reader -> GetOutputPort());
	
	if(UNIT_IN_METERS){
		pgm_filter -> SetScaleMeters();
		mass_properties -> SetScaleMeters();
	} else{
		pgm_filter -> SetScaleKiloMeters();
		mass_properties -> SetScaleKiloMeters();
	}
	pgm_filter -> Update();
	mass_properties -> Update();

	// An instance of SBGATPolyhedronGravityModelUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	SBGATPolyhedronGravityModelUQ pgm_uq;
	pgm_uq.SetPGM(pgm_filter);

	// Populate the shape vertices covariance
	pgm_uq.ComputeVerticesCovarianceGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);

	// Save the covariance
	pgm_uq.SaveNonZeroVerticesCovariance(OUTPUT_DIR + "shape_covariance.json");

	// That's where the grid search should start. Loop over each point in the grid
	// First, create the grid from the bounding box
	std::vector<arma::vec::fixed<3> > grid;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	mass_properties -> GetBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	// Inflate
	xmin *= 2;
	xmax *= 2;
	ymin *= 2;
	ymax *= 2;
	zmin *= 2;
	zmax *= 2;

	// Define grid indices
	int i_max,j_max;

	if (PROJECTION_AXIS == 0){
		i_max = 1./STEP_SIZE * (ymax - ymin);
		j_max = 1./STEP_SIZE * (zmax - zmin);

	}
	else if (PROJECTION_AXIS == 1){

		i_max = 1./STEP_SIZE * (zmax - zmin);
		j_max = 1./STEP_SIZE * (xmax - xmin);

	}
	else if (PROJECTION_AXIS == 2){

		i_max = 1./STEP_SIZE * (xmax - xmin);
		j_max = 1./STEP_SIZE * (zmax - zmin);

	}
	else{
		throw(std::runtime_error("PROJECTION_AXIS has to be 0,1 or 2. It can't be equal to " + std::to_string(PROJECTION_AXIS)));
	}

	// Populate the grid with points that are strictly outside of the reference shape
	for (int i = 0; i < i_max; ++i){
		for (int j = 0; j < j_max; ++j){

			arma::vec::fixed<3> point;
			if (PROJECTION_AXIS == 0){

				double y = ymin + i * STEP_SIZE;
				double z = zmin + j * STEP_SIZE;


				point(0) = 0;
				point(1) = y;
				point(2) = z;
			}
			else if (PROJECTION_AXIS == 1){


				double x = xmin + i * STEP_SIZE;
				double z = zmin + j * STEP_SIZE;

				point(0) = z;
				point(1) = 0;
				point(2) = x;
			}
			else{

				double x = xmin + i * STEP_SIZE;
				double y = ymin + j * STEP_SIZE;

				point(0) = x;
				point(1) = y;
				point(2) = 0;
			}
			
			if (!pgm_filter -> Contains(point))
				grid.push_back(point);
		}
	}
	std::cout << "- Grid size: " << grid.size() << std::endl;

	arma::mat grid_arma(3,grid.size());
	for (int p = 0; p < grid.size(); ++p){
		grid_arma.col(p) = grid[p];
	}
	grid_arma.save(OUTPUT_DIR + "grid.txt",arma::raw_ascii);

	std::cout << "- Sampling grid ...\n";

	// Evaluate the uncertainty at each point on the grid
	arma::vec trace_sqrt_cov_vector(grid.size());

	auto start = std::chrono::system_clock::now();
	boost::progress_display progress(grid.size());
	#pragma omp parallel for
	for (int p = 0; p < grid.size(); ++p){
		trace_sqrt_cov_vector(p) = std::sqrt(arma::trace(pgm_uq.GetCovarianceAcceleration(grid[p])));
		++progress;
	}

	auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
	std::cout << "\n-- Done sampling grid " << elapsed_seconds.count() << " s\n";

	trace_sqrt_cov_vector.save(OUTPUT_DIR + "trace_sqrt_cov_vector.txt",arma::raw_ascii);


	return 0;
}
