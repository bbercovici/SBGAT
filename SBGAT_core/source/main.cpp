#include <iostream>
#include <armadillo>
#include <chrono>

#include "ShapeModelImporter.hpp"
#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"
#include "Constants.hpp"
#include "SyntheticObservations.hpp"

int main( int argc, char** argv ) {

	// Ref frame graph
	SBGAT_CORE::FrameGraph frame_graph;
	frame_graph.add_frame("N");
	frame_graph.add_frame("T");
	frame_graph.add_frame("E");

	frame_graph.add_transform("N", "T");
	frame_graph.add_transform("N", "E");

	SBGAT_CORE::ShapeModel shape_model("T", &frame_graph);
	SBGAT_CORE::ShapeModelImporter shape_io("../../resources/bennu.obj", 1000);

	shape_io.load_shape_model(&shape_model);

	SBGAT_CORE::DynamicAnalyses dynamic_analyses(&shape_model);

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();

	double p[3];
	p[0] = 0;
	p[1] = 0;
	p[2] = 2000;

	std::cout << dynamic_analyses.pgm_acceleration(p, SBGAT_CORE::constants::density::kw4_alpha).t() << std::endl;
	
	
	end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	return 0;
}
