/**
@file   main.cpp
@Author Benjamin Bercovici (bebe0705@colorado.edu)
@date   July, 2017
@brief main file of FrameAccelerations example
*/


#include <iostream>
#include <armadillo>
#include <chrono>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>
#include <Constants.hpp>
#include <RigidBodyKinematics.hpp>

int main( int argc, char** argv ) {

	// A Reference frame graph is created. This is the
	// keystone of SBGAT as it connects different reference frame together
	// and enables coordinates transformations from one frame to another
	SBGAT_CORE::FrameGraph frame_graph;
	frame_graph.add_frame("original");
	frame_graph.add_frame("principal_barycentric");


	// The ECI frame is defined as the parent of the principal_barycentric frame
	// The BF frame and the ECI frame are considered as coincident by default
	frame_graph.add_transform("original", "principal_barycentric");

	// An empty shape model is created
	// Its name is left empty for the moment
	SBGAT_CORE::ShapeModel shape_model("", &frame_graph);

	// A shape model importer is created to load in an OBJ shape model
	SBGAT_CORE::ShapeModelImporter shape_io("../cube.obj", 1);

	// The shape model is built by connecting the shape model importer and the
	// shape model
	shape_io.load_shape_model(&shape_model);

	// The conversion pipeline between frames is created
	arma::vec mrp = RBK::dcm_to_mrp(shape_model.get_original_to_principal_dcm());
	arma::vec cm = shape_model.get_center_of_mass();

	frame_graph.set_transform_mrp("original", "principal_barycentric", mrp);
	frame_graph.set_transform_origin("original", "principal_barycentric", cm);


	// Dynamic analyses are performed on the shape model
	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	// The acceleration of gravity is evaluated at the point of coordinates (1,1,1)
	// in the principal_barycentric frame.
	arma::vec p = {1, 1, 1};

	// The acceleration at the provided point is calculated using the polyhedron gravity model
	arma::vec acc = dynamic_analyses.pgm_acceleration(p.colptr(0), 5000);

	// The acceleration of gravity is evaluated at the point of coordinates (1.5,1.5,1.5)
	// in original frame (that is, the frame that was used in the obj file storing the shape's coordinates)
	arma::vec p_original = {1.5, 1.5, 1.5};

	// The coordinates of p_original are converted to the barycentric frame
	arma::vec p_original_converted = frame_graph.convert(p_original, "original", "principal_barycentric");

	// The acceleration at the provided point is calculated using the polyhedron gravity model
	arma::vec acc_original_converted = dynamic_analyses.pgm_acceleration(p_original_converted.colptr(0), 5000);


	std::cout << "Coordinates of point already in barycentric principal frame: " << std::endl;
	std::cout << p.t() << std::endl;


	std::cout << "Transformed coordinates from original obj-file frame to barycentric principal: " << std::endl;
	std::cout << p_original_converted.t() << std::endl;


	std::cout << "Acceleration from point already in principal barycentric frame: " << std::endl;
	std::cout << acc.t() << std::endl;


	std::cout << "Acceleration from point converted to in principal barycentric frame: " << std::endl;
	std::cout << acc_original_converted.t() << std::endl;

	if (arma::norm(p - p_original_converted) / arma::norm(p_original_converted) < 1e-6) {
		std::cout << "p and p_original_converted are identical!" << std::endl;
	}

	if (arma::norm(acc - acc_original_converted) / arma::norm(acc_original_converted) < 1e-6) {
		std::cout << "The accelerations evaluated at p and p_original_converted are identical!" << std::endl;
	}


	return 0;
}
