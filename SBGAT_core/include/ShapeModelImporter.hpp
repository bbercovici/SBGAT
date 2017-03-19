
#ifndef HEADER_SHAPEMODELIMPORTER
#define HEADER_SHAPEMODELIMPORTER

#include "ShapeModel.hpp"
#include <armadillo>
#include "omp.h"
#include <boost/progress.hpp>


class ShapeModelImporter{

public:
	ShapeModelImporter(std::string filename);
	void load_shape_model(ShapeModel * shape_model, bool construct_edges = true) const;

protected:
	std::string filename;

};

#endif