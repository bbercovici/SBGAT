#include "VoxelGrid.hpp"

VoxelGrid::VoxelGrid(ShapeModel * shape_model, unsigned int N) {


	// The construction of the voxel grid begins by partitioning the space possibly
	// occupied by the shape by a regular 3D grid of Voxel nodes
	// The starting point is the shape's bounding box
	double bounding_box[6];
	shape_model -> get_bounding_box(bounding_box);

	// The bounding box is slightly inflated to make sure that all the vertices are inside
	// the 3D space partition
	bounding_box[0] *= 1.1;
	bounding_box[1] *= 1.1;
	bounding_box[2] *= 1.1;
	bounding_box[3] *= 1.1;
	bounding_box[4] *= 1.1;
	bounding_box[5] *= 1.1;


	// The space is then filled with N * N * N voxel nodes, spanning the
	// inflated bounding box
	std::vector < std::vector<std::vector<std::shared_ptr<VoxelNode> > > > nodes;


	for (unsigned int x_index = 0; x_index < N; ++x_index) {

		std::vector < std::vector<std::shared_ptr<VoxelNode> > > nodes_x_row;

		for (unsigned int y_index = 0; y_index < N; ++y_index) {

			std::vector<std::shared_ptr<VoxelNode> >  nodes_y_row;

			for (unsigned int z_index = 0; z_index < N; ++z_index) {

				double coordinates[3];
				coordinates[0] = bounding_box[0] + x_index * (bounding_box[3] - bounding_box[0]) / (N - 1);
				coordinates[1] = bounding_box[1] + y_index * (bounding_box[4] - bounding_box[1]) / (N - 1);
				coordinates[2] = bounding_box[2] + z_index * (bounding_box[5] - bounding_box[2]) / (N - 1);


				std::shared_ptr<VoxelNode> node = std::make_shared<VoxelNode>(VoxelNode(coordinates));
				nodes_y_row.push_back(node);
			}

			nodes_x_row.push_back(nodes_y_row);
		}

		nodes.push_back(nodes_x_row);

	}

	// Now that this sparse partitioning of the 3D space is available,
	// Voxels are similarly constructed.
	// Note that there are only (N - 1) * (N - 1) * (N - 1) voxels

	std::vector < std::vector<std::vector<std::shared_ptr<Voxel> > > > voxels;

	for (unsigned int x_index = 0; x_index < N - 1; ++x_index) {

		std::vector < std::vector<std::shared_ptr<Voxel> > > voxels_x_row;

		for (unsigned int y_index = 0; y_index < N - 1; ++y_index) {

			std::vector<std::shared_ptr<Voxel> >  voxels_y_row;

			for (unsigned int z_index = 0; z_index < N - 1; ++z_index) {

				std::map<std::string, std::shared_ptr<VoxelNode> > voxel_nodes;

				voxel_nodes["node_0"] = nodes[x_index + 1][y_index][z_index];

				voxel_nodes["node_1"] = nodes[x_index + 1][y_index + 1][z_index];

				voxel_nodes["node_2"] = nodes[x_index][y_index + 1][z_index];

				voxel_nodes["node_3"] = nodes[x_index][y_index][z_index];

				voxel_nodes["node_4"] = nodes[x_index + 1][y_index][z_index + 1];

				voxel_nodes["node_5"] = nodes[x_index + 1][y_index + 1][z_index + 1];

				voxel_nodes["node_6"] = nodes[x_index][y_index + 1][z_index + 1];

				voxel_nodes["node_7"] = nodes[x_index][y_index][z_index + 1];

				std::shared_ptr<Voxel> voxel = std::make_shared<Voxel>(Voxel(voxel_nodes));
				voxels_y_row.push_back(voxel);
			}

			voxels_x_row.push_back(voxels_y_row);
		}

		voxels.push_back(voxels_x_row);

	}

	// Finally, all the voxels are tested for potential inclusion. If a vertex is found 
	// to be empty, it is not added to the voxel grid. This will free the 
	// memory allocated to this voxel. The same thing will happen
	// For the voxel nodes that were not in use



















}