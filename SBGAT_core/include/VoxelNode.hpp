#ifndef HEADER_VOXELNODE
#define HEADER_VOXELNODE

#include "Voxel.hpp"

class Voxel;


class VoxelNode{

public:
	VoxelNode(double * coordinates);
	VoxelNode(Voxel * connected_voxel);

protected:
	std::vector<Voxel *> connected_voxels;
	std::shared_ptr<arma::vec> coordinates;


};

#endif