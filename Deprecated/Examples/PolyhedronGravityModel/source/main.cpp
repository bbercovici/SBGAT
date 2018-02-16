#include <iostream>
#include <armadillo>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>

int main( int argc, char** argv ) {

	// A Reference frame graph is created. This is the
	// keystone of SBGAT as it connects different reference frame together
	// and enables coordinates transformations from one frame to another
	SBGAT_CORE::FrameGraph frame_graph;
	frame_graph.add_frame("ECI");

	// An empty shape model associated with the BF frame is created
	// Note that the frame graph will not be used at all because
	// all operations are performed in the principal body-frame
	SBGAT_CORE::ShapeModel shape_model("BF", &frame_graph);


	// A shape model importer is created to load in an OBJ shape model
	// The vertices coordinates in this OBJ file will be inflated by 1000
	// in order to have the shape internally represented in meters
	SBGAT_CORE::ShapeModelImporter shape_io("../KW4Alpha.obj", 1000);

	// The shape model is built by connecting the shape model importer and the
	// shape model
	// At this stage, the coordinates of the vertices, facet normals, facet centers,...
	// are all expressed in the shape's barycentric principal frame
	shape_io.load_shape_model(&shape_model);

	// Dynamic analyses are performed on the shape model
	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	// First of all, the Polyhedron Gravity Model accelerations are obtained at the center of each facet
	// The resulting accelerations will be saved within each facet of the shape model in a
	// FacetResults object
	double density = 2000; //(kg/m^3)
	double mu = arma::datum::G * density * shape_model . get_volume();

	std::chrono::time_point<std::chrono::system_clock> start, end;
	std::chrono::duration<double> elapsed_seconds ;
	start = std::chrono::system_clock::now();

	dynamic_analyses.compute_pgm_accelerations(mu);

	
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	
	std::cout << "Elapsed time in PGM computation: " << elapsed_seconds.count() << " s" << std::endl;




	return 0;
}
