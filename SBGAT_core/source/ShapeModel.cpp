#include "ShapeModel.hpp"
#include <chrono>

// void ShapeModel::load(const std::string & filename) {

// 	Assimp::Importer importer;

// 	const aiScene * scene = importer.ReadFile( filename, aiProcess_Triangulate | aiProcess_JoinIdenticalVertices  );

// 	// If the import succeeded:
// 	if (scene != NULL) {

// 		// If the imported shape model had at least one mesh in it
// 		if (scene -> mMeshes > 0) {

// 			// For now, only the first mesh is used
// 			this -> vertices = arma::mat(3, scene -> mMeshes[0] -> mNumVertices);
// 			this -> facet_vertices = arma::umat(3, scene -> mMeshes[0] -> mNumFaces);
// 			this -> facet_normals = arma::mat(3, scene -> mMeshes[0] -> mNumFaces);
// 			this -> F_dyads = arma::cube(scene -> mMeshes[0] -> mNumFaces, 3, 3);

// 			this -> NFacets = scene -> mMeshes[0] -> mNumFaces;
// 			this -> NVertices = scene -> mMeshes[0] -> mNumVertices;

// 			std::map<unsigned int , std::set<unsigned int> > vertex_index_to_facet;
// 			std::set<std::set<unsigned int> > edges;

// 			// Vertex coordinates



// 			for (unsigned int vertex = 0; vertex < scene -> mMeshes[0] -> mNumVertices; ++vertex) {

// 				arma::vec vertex_coords = {scene -> mMeshes[0] -> mVertices[vertex].x,
// 				                           scene -> mMeshes[0] -> mVertices[vertex].y,
// 				                           scene -> mMeshes[0] -> mVertices[vertex].z
// 				                          };
// 				this -> vertices.col(vertex) = vertex_coords;
// 			}


// 			// Connectivity Table
// 			for (unsigned int facet = 0; facet < this -> NFacets ; ++facet) {

// 				if (scene -> mMeshes[0] -> mFaces[facet].mNumIndices != 3) {
// 					std::cout << " More than three vertices belong to this facet " << std::endl;
// 					throw " More than three vertices belong to this facet ";
// 				}

// 				this -> facet_vertices.col(facet)(0) = scene -> mMeshes[0] -> mFaces[facet].mIndices[0];
// 				this -> facet_vertices.col(facet)(1) = scene -> mMeshes[0] -> mFaces[facet].mIndices[1];
// 				this -> facet_vertices.col(facet)(2) = scene -> mMeshes[0] -> mFaces[facet].mIndices[2];


// 				vertex_index_to_facet[this -> facet_vertices.col(facet)(0)].insert(facet);
// 				vertex_index_to_facet[this -> facet_vertices.col(facet)(1)].insert(facet);
// 				vertex_index_to_facet[this -> facet_vertices.col(facet)(2)].insert(facet);

// 				std::set<unsigned int> edge_0;
// 				edge_0.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[0]);
// 				edge_0.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[1]);

// 				std::set<unsigned int> edge_1;
// 				edge_1.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[0]);
// 				edge_1.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[2]);

// 				std::set<unsigned int> edge_2;
// 				edge_2.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[1]);
// 				edge_2.insert(scene -> mMeshes[0] -> mFaces[facet].mIndices[2]);

// 				if (this -> edges_to_facets.find(edge_0) == this -> edges_to_facets.end()) {
// 					this -> edges_to_facets[edge_0].insert(facet);
// 				}

// 				else if (this -> edges_to_facets[edge_0].size() < 2) {
// 					this -> edges_to_facets[edge_0].insert(facet);
// 				}

// 				if (this -> edges_to_facets.find(edge_1) == this -> edges_to_facets.end()) {
// 					this -> edges_to_facets[edge_1].insert(facet);
// 				}

// 				else if (this -> edges_to_facets[edge_1].size() < 2) {
// 					this -> edges_to_facets[edge_1].insert(facet);
// 				}


// 				if (this -> edges_to_facets.find(edge_2) == this -> edges_to_facets.end()) {
// 					this -> edges_to_facets[edge_2].insert(facet);
// 				}

// 				else if (this -> edges_to_facets[edge_2].size() < 2) {
// 					this -> edges_to_facets[edge_2].insert(facet);
// 				}

// 				edges.insert(edge_0);
// 				edges.insert(edge_1);
// 				edges.insert(edge_2);

// 			}

// 			this -> NEdges = edges.size();
// 			this -> E_dyads = arma::cube(this -> NEdges, 3, 3);
// 			unsigned int edge_index = 0;


// 			for (std::set<std::set<unsigned int> >::iterator iter = edges.begin(); iter != edges.end(); ++iter) {
// 				this -> edges_to_edges_index[*iter] = edge_index;
// 				this -> edges_indices_to_edge[edge_index] = *iter;
// 				++edge_index;
// 			}

// 			// Normals
// 			#pragma omp parallel for
// 			for (unsigned int facet = 0; facet < this -> NFacets; ++facet) {
// 				unsigned int P0_index = this -> facet_vertices.col(facet)(0);
// 				unsigned int P1_index = this -> facet_vertices.col(facet)(1);
// 				unsigned int P2_index = this -> facet_vertices.col(facet)(2);

// 				arma::vec P0 = this -> vertices.col(P0_index);
// 				arma::vec P1 = this -> vertices.col(P1_index);
// 				arma::vec P2 = this -> vertices.col(P2_index);
// 				arma::vec facet_normal = arma::cross(P1 - P0, P2 - P0) / arma::norm(arma::cross(P1 - P0, P2 - P0));
// 				this -> facet_normals.col(facet) = facet_normal;
// 			}


// 			this -> check_normals_consistency();

// 			this -> compute_dyads();

// 		}


// 	}

// 	else {
// 		std::cout << " There was an error opening the shape model file " << std::endl;
// 		throw " There was an error opening the shape model file ";
// 	}

// }

void ShapeModel::compute_normals() {

	// Normals
	#pragma omp parallel for
	for (unsigned int facet = 0; facet < this -> NFacets; ++facet) {
		unsigned int P0_index = this -> facet_vertices.col(facet)(0);
		unsigned int P1_index = this -> facet_vertices.col(facet)(1);
		unsigned int P2_index = this -> facet_vertices.col(facet)(2);
		arma::vec P0 = this -> vertices.col(P0_index);
		arma::vec P1 = this -> vertices.col(P1_index);
		arma::vec P2 = this -> vertices.col(P2_index);
		arma::vec facet_normal = arma::cross(P1 - P0, P2 - P0) / arma::norm(arma::cross(P1 - P0, P2 - P0));

		this -> facet_normals.col(facet) = facet_normal;
	}

	this -> check_normals_consistency();

}

void ShapeModel::construct_edges() {

	std::map<unsigned int , std::set<unsigned int> > vertex_index_to_facet;
	std::set<std::set<unsigned int> > edges;

	for (unsigned int facet = 0; facet < this -> NFacets ; ++facet) {

		vertex_index_to_facet[this -> facet_vertices.col(facet)(0)].insert(facet);
		vertex_index_to_facet[this -> facet_vertices.col(facet)(1)].insert(facet);
		vertex_index_to_facet[this -> facet_vertices.col(facet)(2)].insert(facet);

		std::set<unsigned int> edge_0;
		edge_0.insert(this -> facet_vertices.col(facet)(0));
		edge_0.insert(this -> facet_vertices.col(facet)(1));

		std::set<unsigned int> edge_1;
		edge_1.insert(this -> facet_vertices.col(facet)(0));
		edge_1.insert(this -> facet_vertices.col(facet)(2));

		std::set<unsigned int> edge_2;
		edge_2.insert(this -> facet_vertices.col(facet)(1));
		edge_2.insert(this -> facet_vertices.col(facet)(2));

		if (this -> edges_to_facets.find(edge_0) == this -> edges_to_facets.end()) {
			this -> edges_to_facets[edge_0].insert(facet);
		}

		else if (this -> edges_to_facets[edge_0].size() < 2) {
			this -> edges_to_facets[edge_0].insert(facet);
		}

		if (this -> edges_to_facets.find(edge_1) == this -> edges_to_facets.end()) {
			this -> edges_to_facets[edge_1].insert(facet);
		}

		else if (this -> edges_to_facets[edge_1].size() < 2) {
			this -> edges_to_facets[edge_1].insert(facet);
		}


		if (this -> edges_to_facets.find(edge_2) == this -> edges_to_facets.end()) {
			this -> edges_to_facets[edge_2].insert(facet);
		}

		else if (this -> edges_to_facets[edge_2].size() < 2) {
			this -> edges_to_facets[edge_2].insert(facet);
		}

		edges.insert(edge_0);
		edges.insert(edge_1);
		edges.insert(edge_2);

	}

	this -> NEdges = edges.size();
	this -> E_dyads = arma::cube(this -> NEdges, 3, 3);

	this -> edges_vertices_indices = arma::umat(2, this -> NEdges);

	for (std::set<std::set<unsigned int> >::iterator iter = edges.begin(); iter != edges.end(); ++iter) {
		this -> edges_to_edges_index[*iter] = std::distance(edges.begin(), iter);
		this -> edges_indices_to_edge.push_back(*iter);
		arma::uvec edge_vertices_indices = {*(iter -> begin()), *std::next(iter -> begin())};
		this -> edges_vertices_indices.col(std::distance(edges.begin(), iter)) = edge_vertices_indices;
	}

}


std::set<unsigned int> ShapeModel::get_edge_from_edge_index(unsigned int edge_index) const {
	return this -> edges_indices_to_edge[edge_index];
}


unsigned int ShapeModel::get_vertex_global_index_from_edge_index(unsigned int vertex_local_index,unsigned int edge_index) const {
	return this -> edges_vertices_indices.col(edge_index)(vertex_local_index);
}


unsigned int ShapeModel::get_NFacets() const {
	return this -> NFacets;
}

unsigned int ShapeModel::get_NVertices() const {
	return this -> NVertices;
}

unsigned int ShapeModel::get_NEdges() const {
	return this -> NEdges;
}

void ShapeModel::compute_F_dyads() {

	for (unsigned int facet = 0; facet < this -> NFacets; ++facet) {
		this -> compute_F_dyad(facet);
	}

}

void ShapeModel::compute_E_dyads() {

	for (auto const & edge : this -> edges_to_facets) {

		this -> compute_E_dyad(edge);

	}

}

void ShapeModel::compute_E_dyad(const std::pair < std::set<unsigned int> ,
                                std::set< unsigned int > > & edge) {

	unsigned int facet_A_index = *edge.second.begin();
	unsigned int facet_B_index = *std::next(edge.second.begin());

	arma::vec facet_normal_A = this -> facet_normals.col(facet_A_index);
	arma::vec facet_normal_B = this -> facet_normals.col(facet_B_index);

	arma::vec edge_direction_vector = arma::normalise(arma::cross(facet_normal_A, facet_normal_B));
	arma::vec edge_normal_A_to_B = arma::normalise(arma::cross(edge_direction_vector, facet_normal_A));
	arma::vec edge_normal_B_to_A = arma::normalise(arma::cross(edge_direction_vector, facet_normal_B));

	// Check if edge normals are consistently oriented
	// In order to do this, the vertices lying on the edge are found
	// along with those that do not belong to the edge

	std::set<unsigned int> vertices_on_edge;
	unsigned int vertex_in_A_global_not_on_edge = 0;
	unsigned int vertex_in_B_global_not_on_edge = 0;

	for (unsigned int vertex_in_A_local = 0; vertex_in_A_local < 3; ++vertex_in_A_local) {

		unsigned int vertex_in_A_global = this -> facet_vertices.col(facet_A_index)(vertex_in_A_local);

		for (unsigned int vertex_in_B_local = 0; vertex_in_B_local < 3; ++vertex_in_B_local) {
			unsigned int vertex_in_B_global = this -> facet_vertices.col(facet_B_index)(vertex_in_B_local);

			if (vertex_in_A_global == vertex_in_B_global) {
				vertices_on_edge.insert(vertex_in_A_global);
			}

		}

	}

	for (unsigned int vertex_in_A_local = 0; vertex_in_A_local < 3; ++vertex_in_A_local) {
		unsigned int vertex_in_A_global = this -> facet_vertices.col(facet_A_index)(vertex_in_A_local);
		if (vertices_on_edge.find(vertex_in_A_global) == vertices_on_edge.end()) {
			vertex_in_A_global_not_on_edge = vertex_in_A_global;
			break;
		}
	}


	for (unsigned int vertex_in_B_local = 0; vertex_in_B_local < 3; ++vertex_in_B_local) {
		unsigned int vertex_in_B_global = this -> facet_vertices.col(facet_B_index)(vertex_in_B_local);
		if (vertices_on_edge.find(vertex_in_B_global) == vertices_on_edge.end()) {
			vertex_in_B_global_not_on_edge = vertex_in_B_global;
			break;
		}
	}

	if (arma::dot(edge_normal_A_to_B,
	              this -> vertices.col(*vertices_on_edge.begin()) - this -> vertices.col(vertex_in_A_global_not_on_edge)) < 0) {
		edge_normal_A_to_B = - edge_normal_A_to_B;
	}


	if (arma::dot(edge_normal_B_to_A,
	              this -> vertices.col(*vertices_on_edge.begin()) - this -> vertices.col(vertex_in_B_global_not_on_edge)) < 0) {
		edge_normal_B_to_A = - edge_normal_B_to_A;
	}


	this -> E_dyads(
	    arma::span(this -> edges_to_edges_index[edge.first]),
	    arma::span(),
	    arma::span()) = facet_normal_A * edge_normal_A_to_B.t() + facet_normal_B * edge_normal_B_to_A.t();

}

void ShapeModel::check_normals_consistency(double tol) const {
	double facet_area_average = 0;

	double sx = 0;
	double sy = 0;
	double sz = 0;

	#pragma omp parallel for reduction(+:facet_area_average,sx,sy,sz)
	for (unsigned int facet = 0; facet < this -> NFacets; ++facet) {

		unsigned int P0_index = this -> facet_vertices.col(facet)(0);
		unsigned int P1_index = this -> facet_vertices.col(facet)(1);
		unsigned int P2_index = this -> facet_vertices.col(facet)(2);

		arma::vec P0 = this -> vertices.col(P0_index);
		arma::vec P1 = this -> vertices.col(P1_index);
		arma::vec P2 = this -> vertices.col(P2_index);
		double facet_area = arma::norm( arma::cross(P1 - P0, P2 - P0)) / 2;

		sx += facet_area * this -> facet_normals.col(facet)(0);
		sy += facet_area * this -> facet_normals.col(facet)(1);
		sz += facet_area * this -> facet_normals.col(facet)(2);
		facet_area_average += facet_area;
	}

	arma::vec surface_sum = {sx, sy, sz};

	facet_area_average = facet_area_average / this -> NFacets;

	if (arma::norm(surface_sum) / facet_area_average > tol) {
		std::cout << "Sum of oriented normals: " << arma::norm(surface_sum) / facet_area_average << std::endl;
		throw "Normals were incorrectly oriented";
	}

}

void ShapeModel::set_NFacets(unsigned int nfacet)  {
	this -> NFacets = nfacet;
}


void ShapeModel::set_NVertices(unsigned int nvertices)  {
	this -> NVertices = nvertices;
}


void ShapeModel::set_NEdges(unsigned int nedges)  {
	this -> NEdges = nedges;
}


void ShapeModel::compute_F_dyad(unsigned int facet) {

	arma::vec normal = this -> facet_normals.col(facet);
	arma::mat facet_dyad = normal * normal.t();

	this -> F_dyads(
	    arma::span(facet),
	    arma::span(),
	    arma::span()) = facet_dyad;
}

void ShapeModel::compute_dyads() {

	this -> compute_F_dyads();
	this -> compute_E_dyads();

}

arma::mat ShapeModel::get_F_dyad(unsigned int facet) const {
	return this -> F_dyads(
	           arma::span(facet),
	           arma::span(),
	           arma::span());
}


void ShapeModel::set_vertices(arma::mat vertices) {
	this -> vertices = vertices;
};

void ShapeModel::set_facet_vertices(arma::umat facet_vertices) {
	this -> facet_vertices = facet_vertices;
	this -> facet_normals = arma::mat(3, this -> NFacets);
	this -> F_dyads = arma::cube(this -> NFacets, 3, 3);

}

void ShapeModel::set_F_dyad(unsigned int facet, arma::mat dyad) {
	this -> F_dyads(
	    arma::span(facet),
	    arma::span(),
	    arma::span()) = dyad;
}

arma::mat ShapeModel::get_E_dyad(unsigned int edge_index) const {
	return this -> E_dyads(
	           arma::span(edge_index),
	           arma::span(),
	           arma::span());
}

void ShapeModel::set_E_dyad(unsigned int edge_index, arma::mat & dyad) {
	this -> E_dyads(
	    arma::span(edge_index),
	    arma::span(),
	    arma::span()) = dyad;
}


arma::uvec ShapeModel::get_vertex_indices_in_facet(unsigned int facet) const {
	return this -> facet_vertices.col(facet);
}

arma::vec ShapeModel::get_vertex(unsigned int vertex_index) const {
	return this -> vertices.col(vertex_index);
}


void ShapeModel::save(const std::string & filename) const {

}


void ShapeModel::load_normals(const std::string & filename) {

	this -> facet_normals = arma::mat(this -> NFacets, 3);
	this -> facet_normals.load(filename);

}

void ShapeModel::save_normals(const std::string & filename)  const {

	this -> facet_normals.save(filename);

}

void ShapeModel::load_F_dyads(const std::string & filename) {
	this -> F_dyads = arma::cube(this -> NEdges, 3, 3);
	this -> F_dyads.load(filename);

}

void ShapeModel::save_F_dyads(const std::string & filename) const {
	this -> F_dyads.save(filename);
}

void ShapeModel::load_E_dyads(const std::string & filename) {
	this -> E_dyads = arma::cube(this -> NEdges, 3, 3);
	this -> E_dyads.load(filename);
}

void ShapeModel::save_E_dyads(const std::string & filename) const {
	this -> E_dyads.save(filename);

}
