#include "VoxelNode.hpp"

VoxelNode::VoxelNode(double * coordinates) {

	arma::vec cords = {coordinates[0], coordinates[1], coordinates[2]};
	this -> coordinates = std::make_shared<arma::vec>(cords);

}

VoxelNode::VoxelNode(Voxel * connected_voxel) {

}