#include "../include/ShapeModel.hpp"
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
// 			}m


// 			this -> check_normals_consistency();

// 			this -> compute_dyads();

// 		}

// 	}

// 	else {
// 		std::cout << " There was an error opening the shape model file " << std::endl;
// 		throw " There was an error opening the shape model file ";
// 	}

// }

ShapeModel::ShapeModel() {

}

void ShapeModel::shift(arma::vec x) {

	// The vertices are shifted
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		*this -> vertices[vertex_index] -> get_coordinates() = *this -> vertices[vertex_index] -> get_coordinates() + x;

	}

	// The facet centers are shifted
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0;
	        facet_index < this -> get_NFacets();
	        ++facet_index) {

		*this -> facets[facet_index] -> get_facet_center() =  *this -> facets[facet_index] -> get_facet_center() + x;
	}


}

ShapeModel::ShapeModel(std::string ref_frame_name,
                       FrameGraph * frame_graph) {
	this -> frame_graph = frame_graph;
	this -> ref_frame_name = ref_frame_name;
}



void ShapeModel::add_facet(Facet * facet) {
	this -> facets. push_back(facet);
}


std::string ShapeModel::get_ref_frame_name() const {
	return this -> ref_frame_name;
}


void ShapeModel::add_edge(std::shared_ptr<Edge> edge) {
	this -> edges. push_back(edge);

}

void ShapeModel::add_vertex(std::shared_ptr<Vertex> vertex) {
	this -> vertices.push_back(vertex);
}

ShapeModel::~ShapeModel() {
	for (unsigned int facet_index = 0; facet_index < this -> facets.size(); ++ facet_index) {
		delete(this -> facets[facet_index]);
	}
}



unsigned int ShapeModel::get_NFacets() const {
	return this -> facets . size();
}

unsigned int ShapeModel::get_NVertices() const {
	return this -> vertices . size();
}

unsigned int ShapeModel::get_NEdges() const {
	return this -> edges . size();
}

std::vector<std::shared_ptr< Vertex> > * ShapeModel::get_vertices() {
	return &this -> vertices;
}


std::vector<Facet * > * ShapeModel::get_facets() {
	return &this -> facets;
}

std::vector<std::shared_ptr< Edge> > * ShapeModel::get_edges() {
	return &this -> edges;
}


void ShapeModel::check_normals_consistency(double tol) const {
	double facet_area_average = 0;

	double sx = 0;
	double sy = 0;
	double sz = 0;

	#pragma omp parallel for reduction(+:facet_area_average,sx,sy,sz) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0; facet_index < this -> facets.size(); ++facet_index) {

		Facet * facet = this -> facets[facet_index];

		sx += facet -> get_area() * facet -> get_facet_normal() -> at(0);
		sy += facet -> get_area() * facet -> get_facet_normal() -> at(1);
		sz += facet -> get_area() * facet -> get_facet_normal() -> at(2);

		facet_area_average += facet -> get_area();

	}

	arma::vec surface_sum = {sx, sy, sz};

	facet_area_average = facet_area_average / this -> facets.size();

	if (arma::norm(surface_sum) / facet_area_average > tol) {
		std::cout << "Sum of oriented normals: " << arma::norm(surface_sum) / facet_area_average << std::endl;
		throw "Normals were incorrectly oriented";
	}

}



void ShapeModel::compute_volume() {
	double volume = 0;

	#pragma omp parallel for reduction(+:volume) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0;
	        facet_index < this -> facets.size();
	        ++facet_index) {

		std::vector<std::shared_ptr<Vertex > > * vertices = this -> facets[facet_index] -> get_vertices();

		arma::vec * r0 =  vertices -> at(0) -> get_coordinates();
		arma::vec * r1 =  vertices -> at(1) -> get_coordinates();
		arma::vec * r2 =  vertices -> at(2) -> get_coordinates();
		double dv = 1. / 6. * arma::dot(*r0, arma::cross(*r1 - *r0, *r2 - *r0));
		volume = volume + dv;

	}

	this -> volume = volume;
}



void ShapeModel::update_mass_properties() {
	this -> compute_surface_area();

	this -> compute_volume();

	this -> compute_center_of_mass();
}

void ShapeModel::compute_center_of_mass() {
	double c_x = 0;
	double c_y = 0;
	double c_z = 0;
	double volume = this -> get_volume();

	#pragma omp parallel for reduction(+:c_x,c_y,c_z) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0;
	        facet_index < this -> facets.size();
	        ++facet_index) {


		std::vector<std::shared_ptr<Vertex > > * vertices = this -> facets[facet_index] -> get_vertices();

		arma::vec * r0 =  vertices -> at(0) -> get_coordinates();
		arma::vec * r1 =  vertices -> at(1) -> get_coordinates();
		arma::vec * r2 =  vertices -> at(2) -> get_coordinates();

		double * r0d =  vertices -> at(0) -> get_coordinates() -> colptr(0);
		double * r1d =  vertices -> at(1) -> get_coordinates() -> colptr(0);
		double * r2d =  vertices -> at(2) -> get_coordinates() -> colptr(0);

		double dv = 1. / 6. * arma::dot(*r1, arma::cross(*r1 - *r0, *r2 - *r0));

		double dr_x = (r0d[0] + r1d[0] + r2d[0]) / 4.;
		double dr_y = (r0d[1] + r1d[1] + r2d[1]) / 4.;
		double dr_z = (r0d[2] + r1d[2] + r2d[2]) / 4.;

		c_x = c_x + dv * dr_x / volume;
		c_y = c_y + dv * dr_y / volume;
		c_z = c_z + dv * dr_z / volume;

	}

	arma::vec cm = {c_x, c_y, c_z};

	this -> cm =  cm ;

}



double ShapeModel::get_volume() const {
	return this -> volume;
}


double ShapeModel::get_surface_area() const {
	return this -> surface_area;
}


arma::vec * ShapeModel::get_center_of_mass() {
	return &(this -> cm);
}


void ShapeModel::compute_surface_area() {
	double surface_area = 0;

	#pragma omp parallel for reduction(+:surface_area) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0; facet_index < this -> facets.size(); ++facet_index) {

		Facet * facet = this -> facets[facet_index];

		surface_area += facet -> get_area();


	}

	this -> surface_area = surface_area;

}



void ShapeModel::get_bounding_box(double * bounding_box) const {

	double xmin = std::numeric_limits<double>::infinity();
	double ymin = std::numeric_limits<double>::infinity();
	double zmin = std::numeric_limits<double>::infinity();

	double xmax =  - std::numeric_limits<double>::infinity();
	double ymax =  - std::numeric_limits<double>::infinity();
	double zmax =  - std::numeric_limits<double>::infinity();

	#pragma omp parallel for reduction(max : xmax,ymax,zmax),reduction(min : xmin,ymin,zmin)
	for ( unsigned int vertex_index = 0; vertex_index < this -> get_NVertices(); ++ vertex_index) {

		double * vertex_cords = this -> vertices[vertex_index] -> get_coordinates() -> colptr(0);

		if (vertex_cords[0] >= xmax) {
			xmax = vertex_cords[0];
		}
		else if (vertex_cords[0] <= xmin) {
			xmin = vertex_cords[0];
		}

		if (vertex_cords[1] >= ymax) {
			ymax = vertex_cords[1];
		}
		else if (vertex_cords[1] <= ymin) {
			ymin = vertex_cords[1];
		}

		if (vertex_cords[2] >= zmax) {
			zmax = vertex_cords[2];
		}
		else if (vertex_cords[2] <= zmin) {
			zmin = vertex_cords[2];
		}

	}

	bounding_box[0] = xmin;
	bounding_box[1] = ymin;
	bounding_box[2] = zmin;
	bounding_box[3] = xmax;
	bounding_box[4] = ymax;
	bounding_box[5] = zmax;


}




void ShapeModel::set_ref_frame_name(std::string ref_frame_name) {

	this -> ref_frame_name = ref_frame_name;
}

