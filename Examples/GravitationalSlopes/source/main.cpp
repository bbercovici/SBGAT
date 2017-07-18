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
	dynamic_analyses.compute_pgm_accelerations(density);

	// The slopes can now be computed, assuming a rotation period of 2.7645 hours and a spin axis oriented along
	// the +z axis in the body-fixed frame (axis of maximum inertia)
	arma::vec spin_axis = {0, 0, 1};
	double spin_rate = 2 * arma::datum::pi / (3600 * 2.7645);
	arma::vec stats_slopes = dynamic_analyses.compute_gravity_slopes(spin_axis, spin_rate);

	// Some statistics about the slopes are available
	std::cout << "Mean slope: "<< stats_slopes(1) << " deg" << std::endl;
	std::cout << "Min slope: " << stats_slopes(0) << " deg" << std::endl;
	std::cout << "Max slope: " << stats_slopes(2) << " deg" << std::endl;

	// One can also loop over all the facets and get their slopes. This is demonstrated with 
	// the first facet of the shape model
	std::cout << "Slope at the center of the first facet of the shape model: " << shape_model.get_facets() -> at(0) -> get_facet_results() -> get_grav_slope() << " deg\n";





	return 0;
}
