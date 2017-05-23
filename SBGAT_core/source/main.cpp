#include <iostream>
#include "VoxelGrid.hpp"
#include "ShapeModelImporter.hpp"
#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"
#include "Densities.hpp"
#include <chrono>

#include <armadillo>



int main( int argc, char** argv ) {

	// Ref frame graph
	FrameGraph frame_graph;
	frame_graph.add_frame("N");
	frame_graph.add_frame("T");
	frame_graph.add_transform("N", "T");


	ShapeModel shape_model("T", &frame_graph);

	// ShapeModelImporter shape_io("../../resources/KW4Alpha.obj",1);
	ShapeModelImporter shape_io("../../resources/sphere.obj", 1);
	// ShapeModelImporter shape_io("../../resources/pumpkin.obj",1000);

	shape_io.load_shape_model(&shape_model);

	DynamicAnalyses dynamic_analyses(&shape_model);

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	dynamic_analyses.compute_pgm(Density::KW4_ALPHA_DENSITY);
	
	end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	return 0;
}
