#ifndef HEADER_VOXEL
#define HEADER_VOXEL

#include "Vertex.hpp"
#include "Facet.hpp"
#include "VoxelNode.hpp"

#include <map>

class Facet;
class Vertex;
class VoxelNode;


class Voxel {

public:

	/**
	Constructor
	*/
	Voxel(std::map<std::string, std::shared_ptr<VoxelNode> > & node_map);


protected:
	std::shared_ptr<std::vector<Facet * > > facets_inside;
	std::shared_ptr<std::vector<Vertex * > > vertices_inside;

	// Voxel * neighbor_Xp = nullptr;
	// Voxel * neighbor_Xm = nullptr;
	// Voxel * neighbor_Yp = nullptr;
	// Voxel * neighbor_Ym = nullptr;
	// Voxel * neighbor_Zp = nullptr;
	// Voxel * neighbor_Zm = nullptr;


	/**
	Numbering of the nodes inside the voxel

			^Z
			|
			|
			|
		7---|------6
	   /|   |     /|
	  / |   |    / |
	 /  |       /  |
	4----------5   |       --------> Y
	|	3------|-- 2
	|  /	/  |  /
	| /	   /   | /
	|/	  /    |/
	0 ---/---- 1
	    /
	   /
	  X

	*/



	std::shared_ptr<VoxelNode> node_0;
	std::shared_ptr<VoxelNode> node_1;
	std::shared_ptr<VoxelNode> node_2;
	std::shared_ptr<VoxelNode> node_3;
	std::shared_ptr<VoxelNode> node_4;
	std::shared_ptr<VoxelNode> node_5;
	std::shared_ptr<VoxelNode> node_6;
	std::shared_ptr<VoxelNode> node_7;


};


#endif