/**
@file   ShapeModelImporter.txt
@Author Benjamin Bercovici (bebe0705@colorado.edu)
@date   May, 2017
@brief  Declaration of the ShapeModelImporter.txt class holding the method and members required to 
load in a Wavefront file storing a shape model.
*/

#ifndef HEADER_SHAPEMODELIMPORTER
#define HEADER_SHAPEMODELIMPORTER

#include "ShapeModel.hpp"
#include <armadillo>
#include "omp.h"

namespace SBGAT_CORE {


	/**
	Declaration of the ShapeModelImporter class. Holds the method and members required to 
	load in a Wavefront file storing a shape model.
	*/
class ShapeModelImporter {

public:

	/**
	Constructor
	@param filename absolute or relative path to the OBJ file to be read
	@param scaling_factor 1 if provided file uses meters, 1000 if km are used, ...
	@param use_edges True if the shape model to be loaded should be augmented with edges
	*/
	ShapeModelImporter(std::string filename, double scaling_factor,
	                   bool use_edges = true);

	/**
	Reads-in an OBJ file storing the shape model info and sets the field of
	$shape_model to the corresponding values
	@param shape_model pointer to the shape model to receive the read data
	@param enforce_principal_axes if true, the loaded shape
	will be rotated to lign up with its principal axes after
	it has been shited to its barycenter
	*/
	void load_shape_model(ShapeModel * shape_model,bool enforce_principal_axes) const;

protected:
	std::string filename;
	double scaling_factor;
	bool use_edges;

};

}

#endif