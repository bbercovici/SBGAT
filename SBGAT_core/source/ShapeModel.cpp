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


void ShapeModel::shift_to_barycenter() {

	arma::vec x = - (*this -> get_center_of_mass());

	// The vertices are shifted
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		*this -> vertices[vertex_index] -> get_coordinates() = *this -> vertices[vertex_index] -> get_coordinates() + x;

	}

	this -> cm = 0 * this -> cm;

}

void ShapeModel::align_with_principal_axes() {

	arma::vec moments;
	arma::mat axes ;
	double l = std::pow(this -> volume, 1. / 3.);
	arma::mat non_dim_I = this -> body_inertia / (l * l * l * l * l);

	double T = arma::trace(non_dim_I) ;
	double Pi = 0.5 * (T * T - arma::trace(non_dim_I * non_dim_I));
	double U = std::sqrt(T * T - 3 * Pi) / 3;
	double Det = arma::det(non_dim_I);

	if (U > 1e-6) {

		double Theta = std::acos( (- 2 * T * T * T +  9 * T * Pi - 27 * Det) / (54 * U * U * U ));

		double A = T / 3 - 2 * U * std::cos(Theta / 3);
		double B = T / 3 - 2 * U * std::cos(Theta / 3 - 2 * arma::datum::pi / 3);
		double C = T / 3 - 2 * U * std::cos(Theta / 3 + 2 * arma::datum::pi / 3);

		moments = {A, B, C};

		arma::mat L0 = non_dim_I - moments(0) * arma::eye<arma::mat>(3, 3);
		arma::mat L1 = non_dim_I - moments(1) * arma::eye<arma::mat>(3, 3);

		L0.row(0) = arma::normalise(L0.row(0));
		L0.row(1) = arma::normalise(L0.row(1));
		L0.row(2) = arma::normalise(L0.row(2));

		L1.row(0) = arma::normalise(L1.row(0));
		L1.row(1) = arma::normalise(L1.row(1));
		L1.row(2) = arma::normalise(L1.row(2));

		arma::mat e0_mat(3, 3);

		e0_mat.row(0) = arma::cross(L0.row(0), L0.row(1));
		e0_mat.row(1) = arma::cross(L0.row(0), L0.row(2));
		e0_mat.row(2) = arma::cross(L0.row(1), L0.row(2));

		arma::vec norms_e0 = {arma::norm(e0_mat.row(0)), arma::norm(e0_mat.row(1)), arma::norm(e0_mat.row(2))};
		double best_e0 = norms_e0.index_max();
		arma::vec e0 = arma::normalise(e0_mat.row(best_e0).t());

		arma::mat e1_mat(3, 3);
		e1_mat.row(0) = arma::cross(L1.row(0), L1.row(1));
		e1_mat.row(1) = arma::cross(L1.row(0), L1.row(2));
		e1_mat.row(2) = arma::cross(L1.row(1), L1.row(2));

		arma::vec norms_e1 = {arma::norm(e1_mat.row(0)), arma::norm(e1_mat.row(1)), arma::norm(e1_mat.row(2))};
		double best_e1 = norms_e1.index_max();
		arma::vec e1 = arma::normalise(e1_mat.row(best_e1).t());

		arma::vec e2 = arma::cross(e0, e1);

		axes = arma::join_rows(e0, arma::join_rows(e1, e2));
	}

	else {
		moments = std::pow(Det, 1. / 3.) * arma::ones<arma::vec>(3);
		axes = arma::eye<arma::mat>(3, 3);
	}

	moments = l * l * l * l * l * moments;



	std::cout << "Non-dimensional inertia: " << std::endl;
	std::cout << non_dim_I << std::endl;
	
	std::cout << "Principal axes: " << std::endl;
	std::cout << axes << std::endl;

	std::cout << "Principal moments: " << std::endl;
	std::cout << moments << std::endl;




	// The vertices are shifted
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		*this -> vertices[vertex_index] -> get_coordinates() = axes.t() * (*this -> vertices[vertex_index] -> get_coordinates());

	}

	this -> body_inertia = arma::diagmat(moments);
	this -> inertia_axes = axes;

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
	this -> compute_inertia();

}

void ShapeModel::update_facets() {

	for (auto & facet : this -> facets) {
		facet -> update();
	}

}


void ShapeModel::update_edges() {

	for (auto & edge : this -> edges) {
		edge -> compute_dyad();
	}

}

bool ShapeModel::contains(double * point, double tol ) {

	double lagrangian = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:lagrangian) if (USE_OMP_DYNAMIC_ANALYSIS)
	for (unsigned int facet_index = 0; facet_index < this  -> get_NFacets(); ++ facet_index) {

		std::vector<std::shared_ptr<Vertex > > * vertices = this  -> get_facets() -> at(facet_index) -> get_vertices();

		const double * r1 =  vertices -> at(0) -> get_coordinates() -> colptr(0);
		const double * r2 =  vertices -> at(1) -> get_coordinates() -> colptr(0);
		const double * r3 =  vertices -> at(2) -> get_coordinates() -> colptr(0);

		double r1m[3];
		double r2m[3];
		double r3m[3];

		r1m[0] = r1[0] - point[0];
		r1m[1] = r1[1] - point[1];
		r1m[2] = r1[2] - point[2];

		r2m[0] = r2[0] - point[0];
		r2m[1] = r2[1] - point[1];
		r2m[2] = r2[2] - point[2];

		r3m[0] = r3[0] - point[0];
		r3m[1] = r3[1] - point[1];
		r3m[2] = r3[2] - point[2];


		double R1 = std::sqrt( r1m[0] * r1m[0]
		                       + r1m[1] * r1m[1]
		                       + r1m[2] * r1m[2]       );

		double R2 = std::sqrt( r2m[0] * r2m[0]
		                       + r2m[1] * r2m[1]
		                       + r2m[2] * r2m[2]      );


		double R3 = std::sqrt( r3m[0] * r3m[0]
		                       + r3m[1] * r3m[1]
		                       + r3m[2] * r3m[2]      );

		double r2_cross_r3_0 = r2m[1] * r3m[2] - r2m[2] * r3m[1];
		double r2_cross_r3_1 = r3m[0] * r2m[2] - r3m[2] * r2m[0];
		double r2_cross_r3_2 = r2m[0] * r3m[1] - r2m[1] * r3m[0];


		double wf = 2 * std::atan2(
		                r1m[0] * r2_cross_r3_0 + r1m[1] * r2_cross_r3_1 + r1m[2] * r2_cross_r3_2,

		                R1 * R2 * R3 + R1 * (r2m[0] * r3m[0] + r2m[1] * r3m[1]  + r2m[2] * r3m[2] )
		                + R2 * (r3m[0] * r1m[0] + r3m[1] * r1m[1] + r3m[2] * r1m[2])
		                + R3 * (r1m[0] * r2m[0] + r1m[1] * r2m[1] + r1m[2] * r2m[2]));



		lagrangian += wf;

	}

	if (std::abs(lagrangian) < tol) {
		return false;
	}
	else {
		return true;
	}



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

void ShapeModel::compute_inertia() {


	double P_xx = 0;
	double P_yy = 0;
	double P_zz = 0;
	double P_xy = 0;
	double P_xz = 0;
	double P_yz = 0;

	#pragma omp parallel for reduction(+:P_xx,P_yy,P_zz,P_xy,P_xz,P_yz) if (USE_OMP_SHAPE_MODEL)
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



		P_xx += dv / 20 * (2 * r0d[0] * r0d[0]
		                   + 2 * r1d[0] * r1d[0]
		                   + 2 * r2d[0] * r2d[0]
		                   + r0d[0] * r1d[0]
		                   + r0d[0] * r1d[0]
		                   + r0d[0] * r2d[0]
		                   + r0d[0] * r2d[0]
		                   + r1d[0] * r2d[0]
		                   + r1d[0] * r2d[0]);


		P_yy += dv / 20 * (2 * r0d[1] * r0d[1]
		                   + 2 * r1d[1] * r1d[1]
		                   + 2 * r2d[1] * r2d[1]
		                   + r0d[1] * r1d[1]
		                   + r0d[1] * r1d[1]
		                   + r0d[1] * r2d[1]
		                   + r0d[1] * r2d[1]
		                   + r1d[1] * r2d[1]
		                   + r1d[1] * r2d[1]);

		P_zz += dv / 20 * (2 * r0d[2] * r0d[2]
		                   + 2 * r1d[2] * r1d[2]
		                   + 2 * r2d[2] * r2d[2]
		                   + r0d[2] * r1d[2]
		                   + r0d[2] * r1d[2]
		                   + r0d[2] * r2d[2]
		                   + r0d[2] * r2d[2]
		                   + r1d[2] * r2d[2]
		                   + r1d[2] * r2d[2]);

		P_xy += dv / 20 * (2 * r0d[0] * r0d[1]
		                   + 2 * r1d[0] * r1d[1]
		                   + 2 * r2d[0] * r2d[1]
		                   + r0d[0] * r1d[1]
		                   + r0d[1] * r1d[0]
		                   + r0d[0] * r2d[1]
		                   + r0d[1] * r2d[0]
		                   + r1d[0] * r2d[1]
		                   + r1d[1] * r2d[0]);

		P_xz += dv / 20 * (2 * r0d[0] * r0d[2]
		                   + 2 * r1d[0] * r1d[2]
		                   + 2 * r2d[0] * r2d[2]
		                   + r0d[0] * r1d[2]
		                   + r0d[2] * r1d[0]
		                   + r0d[0] * r2d[2]
		                   + r0d[2] * r2d[0]
		                   + r1d[0] * r2d[2]
		                   + r1d[2] * r2d[0]);

		P_yz += dv / 20 * (2 * r0d[1] * r0d[2]
		                   + 2 * r1d[1] * r1d[2]
		                   + 2 * r2d[1] * r2d[2]
		                   + r0d[1] * r1d[2]
		                   + r0d[2] * r1d[1]
		                   + r0d[1] * r2d[2]
		                   + r0d[2] * r2d[1]
		                   + r1d[1] * r2d[2]
		                   + r1d[2] * r2d[1]);

	}

	// The inertia tensor is finally assembled

	arma::mat I = {
		{P_yy + P_zz, -P_xy, -P_xz},
		{ -P_xy, P_xx + P_zz, -P_yz},
		{ -P_xz, -P_yz, P_xx + P_yy}
	};

	this -> body_inertia = I;


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




