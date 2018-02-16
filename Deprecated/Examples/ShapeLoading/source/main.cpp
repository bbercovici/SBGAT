#include <iostream>
#include <armadillo>
#include <chrono>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>
#include <Constants.hpp>

int main( int argc, char** argv ) {

	// A Reference frame graph is created. This is the 
	// keystone of SBGAT as it connects different reference frame together
	// and enables coordinates transformations from one frame to another 
	SBGAT_CORE::FrameGraph frame_graph;
	frame_graph.add_frame("ECI");
	frame_graph.add_frame("BF");
	
	// The ECI frame is defined as the parent of the body-fixed frame
	// The BF frame and the ECI frame are considered as coincident by default
	frame_graph.add_transform("ECI", "BF");

	// An empty shape model associated with the BF frame is created
	SBGAT_CORE::ShapeModel shape_model("BF", &frame_graph);

	// A shape model importer is created to load in an OBJ shape model 
	// The vertices coordinates in this OBJ file will be inflated by 1000
	// in order to have the shape internally represented in meters
	SBGAT_CORE::ShapeModelImporter shape_io("../../GravitationalSlopes/KW4Alpha.obj", 1000);

	// The shape model is built by connecting the shape model importer and the 
	// shape model
	shape_io.load_shape_model(&shape_model);

	// The dimensionless inertia expressed in the body-fixed frame is shown
	std::cout << "Body-fixed inertia tensor (kg m^2):\n";

	std::cout << shape_model.get_inertia() << std::endl;

	// Dynamic analyses are performed on the shape model
	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	// The acceleration of gravity is evaluated at the point of coordinates (0,0,2000)
	// in the body-fixed frame
	double p[3];
	p[0] = 0;
	p[1] = 0;
	p[2] = 2000;

	// The acceleration at the provided point is calculated using the polyhedron gravity model
	double mu = SBGAT_CORE::constants::density::kw4_alpha * arma::datum::G * shape_model . get_volume();
		
	std::cout << "Gravity acceleration at (0,0,2000) (m) in m/s^2:";
	std::cout << dynamic_analyses.pgm_acceleration(p, mu).t() << std::endl;
	


	return 0;
}
