#ifndef HEADER_VOXELGRID
#define HEADER_VOXELGRID

#include "ShapeModel.hpp"
#include "VoxelNode.hpp"
#include "Voxel.hpp"

class Voxel;

class VoxelNode;

class ShapeModel;

class VoxelGrid {

public:

	/**
	Constructor
	@param shape_model Pointer to the shape model that the 
	voxel grid should capture
	@param N Granularity of voxel grid (i.e number of nodes per dimension) 
	*/
	VoxelGrid(ShapeModel * shape_model, unsigned int N = 100);



protected:

	ShapeModel * shape_model;
	// std::vector<std::shared_ptr<Voxel> > voxels;

};

#endif