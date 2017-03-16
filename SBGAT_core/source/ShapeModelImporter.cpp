#include "ShapeModelImporter.hpp"

ShapeModelImporter::ShapeModelImporter(std::string filename) {
	this -> filename = filename;
}

void ShapeModelImporter::load_shape_model(ShapeModel * shape_model, bool construct_edges ) const {

	std::ifstream ifs(this -> filename);

	if (!ifs.is_open()) {
		std::cout << "There was a problem opening the input file!\n";
		throw;
	}

	std::string line;
	std::vector<arma::vec> vertices;
	std::vector<arma::uvec> facet_vertices;

	while (std::getline(ifs, line)) {

		std::stringstream linestream(line);

		std::string word1, word2, word3;

		char type;
		linestream >> type;

		if (type == '#') {
			continue;
		}

		else if (type == 'v') {
			double vx, vy, vz;
			linestream >> vx >> vy >> vz;
			arma::vec vertex = {vx, vy, vz};
			vertices.push_back(vertex);

		}

		else if (type == 'f') {
			unsigned int v0, v1, v2;
			linestream >> v0 >> v1 >> v2;
			arma::uvec vertices_in_facet = {v0 - 1, v1 - 1, v2 - 1};
			facet_vertices.push_back(vertices_in_facet);

		}

		else {
			std::cout << " unrecognized type: " << type << std::endl;
			throw;
		}

	}

	std::cout << " Number of vertices: " << vertices.size() << std::endl;
	std::cout << " Number of facets: " << facet_vertices.size() << std::endl;

	// The std vectors are turned into proper Armadillo types
	arma::mat vertices_arma = arma::mat(3, vertices.size());
	arma::umat facet_vertices_arma = arma::umat(3, facet_vertices.size());

	#pragma omp parallel for
	for (unsigned int vertex = 0; vertex < vertices.size(); ++vertex) {

		arma::vec vertex_coord = {
			vertices[vertex](0),
			vertices[vertex](1),
			vertices[vertex](2)
		};
		vertices_arma.col(vertex) = vertex_coord;

	}

	#pragma omp parallel for
	for (unsigned int facet = 0; facet < facet_vertices.size(); ++facet) {

		arma::uvec facet_vertices_vec = {
			facet_vertices[facet](0),
			facet_vertices[facet](1),
			facet_vertices[facet](2)
		};

		facet_vertices_arma.col(facet) = facet_vertices_vec;

	}

	shape_model -> set_NFacets(facet_vertices.size());
	shape_model -> set_NVertices(vertices.size());
	shape_model -> set_vertices(vertices_arma);
	shape_model -> set_facet_vertices(facet_vertices_arma);

	shape_model -> compute_normals();


	if (construct_edges == true) {

		shape_model -> construct_edges();
		shape_model -> compute_dyads();

	}


}
