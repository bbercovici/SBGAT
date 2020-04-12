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
	bool HOLD_MASS_CONSTANT  = input_data["HOLD_MASS_CONSTANT"];


	int PROJECTION_AXIS = input_data["PROJECTION_AXIS"];
	int N_MONTE_CARLO = input_data["N_MONTE_CARLO"];

	std::string OUTPUT_DIR = input_data["OUTPUT_DIR"];

	std::cout << "- Path to shape: " << PATH_SHAPE << std::endl;
	std::cout << "- Standard deviation on point coordinates (m) : " << ERROR_STANDARD_DEV << std::endl;
	std::cout << "- Correlation distance (m) : " << CORRELATION_DISTANCE << std::endl;
	std::cout << "- Density (kg/m^3) : " << DENSITY << std::endl;
	std::cout << "- Monte-Carlo draws : " << N_MONTE_CARLO << std::endl;
	std::cout << "- Step size : " << STEP_SIZE << std::endl;
	std::cout << "- Projection axis : " << PROJECTION_AXIS << std::endl;


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

	// An instance of SBGATPolyhedronGravityModelUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	SBGATPolyhedronGravityModelUQ pgm_uq;
	pgm_uq.SetModel(pgm_filter);
	pgm_uq.PrecomputeMassPropertiesPartials();
	
	// Populate the shape vertices covariance
	pgm_uq.ComputeVerticesCovarianceGlobal(ERROR_STANDARD_DEV,CORRELATION_DISTANCE);
	// Regularizing the covariance
	int regularized_eigen_values = pgm_uq.RegularizeCovariance();

	std::cout << regularized_eigen_values << " eigenvalues were regularized\n";
	
	arma::mat C_CC = pgm_uq.GetCovarianceSquareRoot();
	arma::mat P_CC = pgm_uq.GetVerticesCovariance();
	std::cout << "Maximum absolute error in covariance square root: " << arma::abs(P_CC - C_CC * C_CC.t()).max() << std::endl;
	
	std::cout << "Saving shape covariance ...\n";
	P_CC.save(OUTPUT_DIR + "full_covariance.txt",arma::raw_ascii);

	// Saving baseline slices
	pgm_uq.TakeAndSaveSlice(0,OUTPUT_DIR + "baseline_slice_x.txt",1e-6);
	pgm_uq.TakeAndSaveSlice(1,OUTPUT_DIR + "baseline_slice_y.txt",1e-6);
	pgm_uq.TakeAndSaveSlice(2,OUTPUT_DIR + "baseline_slice_z.txt",1e-6);

	// First, create the grid from the bounding box
	std::vector<arma::vec::fixed<3> > grid;
	double xmin,xmax,ymin,ymax,zmin,zmax;
	pgm_filter -> GetBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);
	
	// Inflate
	xmin *= 2.5;
	xmax *= 2.5;
	ymin *= 2.5;
	ymax *= 2.5;
	zmin *= 2.5;
	zmax *= 2.5;

	double min_dim = std::min(xmin,std::min(ymin,zmin));
	double max_dim = std::max(xmax,std::max(ymax,zmax));

	xmax = max_dim;
	ymax = max_dim;
	zmax = max_dim;

	xmin = min_dim;
	ymin = min_dim;
	zmin = min_dim;

	// Define grid indices
	int i_max,j_max;

	if (PROJECTION_AXIS == 0){
		i_max = 1./STEP_SIZE * (ymax - ymin);
		j_max = 1./STEP_SIZE * (zmax - zmin);

	}
	else if (PROJECTION_AXIS == 1){

		i_max = 1./STEP_SIZE * (xmax - xmin);
		j_max = 1./STEP_SIZE * (zmax - zmin);

	}
	else if (PROJECTION_AXIS == 2){

		i_max = 1./STEP_SIZE * (xmax - xmin);
		j_max = 1./STEP_SIZE * (ymax - ymin);

	}
	else{
		throw(std::runtime_error("PROJECTION_AXIS has to be either 0, 1 or 2. It can't be equal to " + std::to_string(PROJECTION_AXIS)));
	}

	// Create the grid containers
	std::vector<std::vector<int> > indices;
	arma::mat trace_sqrt_cov(i_max,j_max);
	arma::mat reference_acceleration(i_max,j_max);
	arma::mat uncertainty_over_reference_acc_percentage(i_max,j_max);
	arma::mat inside_outside(i_max,j_max);
	


	// Construct the grid (constant, uniform step size in both directions of space)
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

				point(0) = x;
				point(1) = 0;
				point(2) = z;
			}
			else{

				double x = xmin + i * STEP_SIZE;
				double y = ymin + j * STEP_SIZE;

				point(0) = x;
				point(1) = y;
				point(2) = 0;
			}
			
			grid.push_back(point);
			indices.push_back(std::vector<int>({i,j}));

			if (pgm_filter -> Contains(point)){
				inside_outside(i,j) = 1;
			}
			else{
				inside_outside(i,j) = 0;
			}
		}
	}
	std::cout << "- Grid size: " << grid.size() << std::endl;
	
	// Saving the grid to a file
	arma::mat grid_arma(3,grid.size());
	for (unsigned int p = 0; p < grid.size(); ++p){
		grid_arma.col(p) = grid[p];
	}
	grid_arma.save(OUTPUT_DIR + "grid.txt",arma::raw_ascii);

	// Running a Monte Carlo on a subset of the grid
	std::vector<std::vector<arma::vec::fixed<3> > > all_accelerations;
	std::vector<std::vector<bool>> all_inside_outside;
	std::vector<std::vector<double > > all_potentials;
	std::vector<arma::vec> deviations;
	std::vector<double> densities;


	// arma::vec::fixed<3> e0 = {1,0,0};
	// arma::vec::fixed<3> e1 = {0,1,0};
	// arma::vec::fixed<3> e2 = {0,0,1};
	// arma::vec::fixed<3> e3 = arma::normalise(arma::vec({1,1,0}));
	// arma::vec::fixed<3> e4 = arma::normalise(arma::vec({0,1,1}));
	// arma::vec::fixed<3> e5 = arma::normalise(arma::vec({1,0,1}));
	// arma::vec::fixed<3> e6 = arma::normalise(arma::vec({-1,1,0}));
	// arma::vec::fixed<3> e7 = arma::normalise(arma::vec({0,-1,1}));
	// arma::vec::fixed<3> e8 = arma::normalise(arma::vec({1,0,-1}));


	arma::vec::fixed<3> e0 = {1,0,0};
	arma::vec::fixed<3> e1 = {0,0.8,0};
	arma::vec::fixed<3> e2 = {0,0,0.8};
	arma::vec::fixed<3> e3 = 0.8 * arma::normalise(arma::vec({1,1,0}));
	arma::vec::fixed<3> e4 = 0.8 * arma::normalise(arma::vec({0,1,1}));
	arma::vec::fixed<3> e5 = arma::normalise(arma::vec({1,0,1}));
	arma::vec::fixed<3> e6 = 0.8 * arma::normalise(arma::vec({-1,1,0}));
	arma::vec::fixed<3> e7 = 0.8 * arma::normalise(arma::vec({0,-1,1}));
	arma::vec::fixed<3> e8 = arma::normalise(arma::vec({1,0,-1}));


	std::vector<arma::vec::fixed<3> > all_positions;
	std::vector<double> distances = {200,300,400,500,600};
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



	std::cout << "Running MC ... ";

	auto start = std::chrono::system_clock::now();
	SBGATPolyhedronGravityModelUQ::RunMCUQAccelerationInertial(PATH_SHAPE,DENSITY,
		UNIT_IN_METERS,
		HOLD_MASS_CONSTANT,
		C_CC,
		N_MONTE_CARLO, 
		all_positions,
		OUTPUT_DIR,
		std::min(30,N_MONTE_CARLO),
		deviations,
		densities,
		all_accelerations,
		all_inside_outside);
	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;

	std::cout << "Done running MC in " << elapsed_seconds.count() << " s\n";
		
// Computing MC Dispersions
	arma::vec KL_divergence_analytical_vs_mc(all_positions.size());
	arma::vec abs_value_cov_difference_analytical_vs_mc(all_positions.size());
	arma::vec rel_value_cov_difference_analytical_vs_mc(all_positions.size());
	arma::mat all_positions_arma(3,all_positions.size());
	arma::vec all_inside_outside_arma = arma::zeros<arma::vec>(all_positions.size());


	#pragma omp parallel for
	for (int e = 0; e < all_positions.size(); ++e){

		arma::mat accelerations_mc(3,N_MONTE_CARLO);

		for (int sample = 0; sample < N_MONTE_CARLO; ++sample){
			accelerations_mc.col(sample) = all_accelerations[sample][e];
			if (all_inside_outside[sample][e] == true){
				all_inside_outside_arma(e) += 1;
			}

		}

		arma::mat mc_covariances_acc = arma::cov(accelerations_mc.t());

		arma::vec mc_mean_acc = arma::mean(accelerations_mc,1);
		arma::vec reference_acc = pgm_filter -> GetAcceleration(all_positions[e]);

		arma::mat cov_analytical = pgm_uq.GetCovarianceAcceleration(all_positions[e],HOLD_MASS_CONSTANT);

		KL_divergence_analytical_vs_mc(e) = SBGATFilterUQ::KLDivergence(reference_acc,mc_mean_acc,	cov_analytical,mc_covariances_acc);

		abs_value_cov_difference_analytical_vs_mc(e) = arma::abs(arma::vectorise(cov_analytical - mc_covariances_acc)).max();
		rel_value_cov_difference_analytical_vs_mc(e) = arma::norm(cov_analytical - mc_covariances_acc)/arma::trace(mc_covariances_acc) * 100;
		all_positions_arma.col(e) = all_positions[e];
	}

	


	std::cout << "- Evaluating analytical uncertainties over grid ...\n";

	start = std::chrono::system_clock::now();


	
	// For every point in the grid, evaluate the analytical acceleration covariance and 
	// run a monte carlo to get a statistical covariance to compare against
	boost::progress_display progress(grid.size());

	#pragma omp parallel for
	for (int p = 0; p < grid.size(); ++p){

		int i = indices[p][0];
		int j = indices[p][1];

		const arma::vec::fixed<3> & grid_point = grid[p];

		arma::vec::fixed<3> reference_acceleration_vector = pgm_filter -> GetAcceleration(grid_point);
		arma::mat::fixed<3,3> covariance_acceleration_analytical = pgm_uq.GetCovarianceAcceleration(grid_point,HOLD_MASS_CONSTANT);
		
		reference_acceleration(i,j) = arma::norm(reference_acceleration_vector);
		trace_sqrt_cov(i,j) = std::sqrt(arma::trace(covariance_acceleration_analytical)) ;
		uncertainty_over_reference_acc_percentage(i,j) = trace_sqrt_cov(i,j) / reference_acceleration(i,j) * 100;

		++progress;
	}

	end = std::chrono::system_clock::now();
	elapsed_seconds = end-start;
	std::cout << "\n-- Done evaluating over grid in " << elapsed_seconds.count() << " s\n";

	trace_sqrt_cov.save(OUTPUT_DIR + "trace_sqrt_cov.txt",arma::raw_ascii);
	reference_acceleration.save(OUTPUT_DIR + "reference_acceleration.txt",arma::raw_ascii);
	inside_outside.save(OUTPUT_DIR + "inside_outside.txt",arma::raw_ascii);
	uncertainty_over_reference_acc_percentage.save(OUTPUT_DIR + "uncertainty_over_reference_acc_percentage.txt",arma::raw_ascii);

	abs_value_cov_difference_analytical_vs_mc.save(OUTPUT_DIR + "abs_value_cov_difference_analytical_vs_mc.txt",arma::raw_ascii);
	rel_value_cov_difference_analytical_vs_mc.save(OUTPUT_DIR + "rel_value_cov_difference_analytical_vs_mc.txt",arma::raw_ascii);
	KL_divergence_analytical_vs_mc.save(OUTPUT_DIR + "KL_divergence_analytical_vs_mc.txt",arma::raw_ascii);
	all_positions_arma.save(OUTPUT_DIR + "all_positions_arma.txt",arma::raw_ascii);
	all_inside_outside_arma.save(OUTPUT_DIR + "all_inside_outside_arma.txt",arma::raw_ascii);

	return 0;
}
