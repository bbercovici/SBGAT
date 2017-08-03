#include "ShapeModel.hpp"

using namespace SBGAT_CORE;

ShapeModel::ShapeModel() {

}

void ShapeModel::update_mass_properties() {

	this -> barycenter_aligned = false;
	this -> principal_axes_aligned = false;

	this -> compute_surface_area();
	this -> compute_volume();
	this -> compute_center_of_mass();
	this -> compute_inertia();
	this -> compute_principal_axes();



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

void ShapeModel::add_facet(Facet * facet) {
	this -> facets. push_back(facet);
}



void ShapeModel::update_facets(std::set<Facet *> & facets) {

	for (auto & facet : facets) {
		facet -> update();
	}

}




void ShapeModel::shift_rotate_to_principal_frame() {

	if (this -> barycenter_aligned == false) {
		this -> shift_to_barycenter();
		this -> barycenter_aligned = true;
	}

	if (this -> barycenter_aligned == true &&
	        this -> principal_axes_aligned == false) {
		this -> align_with_principal_axes();
		this -> principal_axes_aligned = true;

	}

}


bool ShapeModel::contains(double * point, double tol ) {

	double laplacian = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:laplacian) if (USE_OMP_DYNAMIC_ANALYSIS)
	for (unsigned int facet_index = 0; facet_index < this -> get_NFacets(); ++ facet_index) {

		std::vector<std::shared_ptr<SBGAT_CORE::Vertex > > * vertices = this -> get_facets() -> at(facet_index) -> get_vertices();

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



		laplacian += wf;

	}

	if (std::abs(laplacian) < tol) {
		return false;
	}
	else {
		return true;
	}



}


void ShapeModel::save(std::string path) const {
	std::ofstream shape_file;
	shape_file.open(path);

	std::map<std::shared_ptr<Vertex> , unsigned int> vertex_ptr_to_index;

	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		shape_file << "v " << this -> vertices[vertex_index] -> get_coordinates() -> colptr(0)[0] << " " << this -> vertices[vertex_index] -> get_coordinates() -> colptr(0)[1] << " " << this -> vertices[vertex_index] -> get_coordinates() -> colptr(0)[2] << std::endl;
		vertex_ptr_to_index[this -> vertices[vertex_index]] = vertex_index;
	}

	for (unsigned int facet_index = 0;
	        facet_index < this -> get_NFacets();
	        ++facet_index) {

		unsigned int v0 =  vertex_ptr_to_index[this -> facets[facet_index] -> get_vertices() -> at(0)] + 1;
		unsigned int v1 =  vertex_ptr_to_index[this -> facets[facet_index] -> get_vertices() -> at(1)] + 1;
		unsigned int v2 =  vertex_ptr_to_index[this -> facets[facet_index] -> get_vertices() -> at(2)] + 1;

		shape_file << "f " << v0 << " " << v1 << " " << v2 << std::endl;

	}




	shape_file.close();




}


void ShapeModel::shift_to_barycenter() {


	arma::vec x = - this -> cm;

	// The vertices are shifted
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		*this -> vertices[vertex_index] -> get_coordinates() = *this -> vertices[vertex_index] -> get_coordinates() + x;

	}


}

void ShapeModel::align_with_principal_axes() {


	// The vertices are rotated so as to have their coordinates
	// expressed in the principal frame
	#pragma omp parallel for if(USE_OMP_SHAPE_MODEL)
	for (unsigned int vertex_index = 0;
	        vertex_index < this -> get_NVertices();
	        ++vertex_index) {

		*this -> vertices[vertex_index] -> get_coordinates() = this -> original_to_principal_dcm * (*this -> vertices[vertex_index] -> get_coordinates());
	}


}

arma::mat ShapeModel::get_inertia() const {
	return this -> inertia;

}



void ShapeModel::compute_principal_axes() {


	arma::vec moments;
	arma::mat axes ;

	double T = arma::trace(this -> inertia) ;
	double Pi = 0.5 * (T * T - arma::trace(this -> inertia * this -> inertia));
	double U = std::sqrt(T * T - 3 * Pi) / 3;
	double Det = arma::det(this -> inertia);

	if (U > 1e-6) {

		double Theta = std::acos( (- 2 * T * T * T +  9 * T * Pi - 27 * Det) / (54 * U * U * U ));

		double A = T / 3 - 2 * U * std::cos(Theta / 3);
		double B = T / 3 - 2 * U * std::cos(Theta / 3 - 2 * arma::datum::pi / 3);
		double C = T / 3 - 2 * U * std::cos(Theta / 3 + 2 * arma::datum::pi / 3);

		moments = {A, B, C};

		arma::mat L0 = this -> inertia - moments(0) * arma::eye<arma::mat>(3, 3);
		arma::mat L1 = this -> inertia - moments(1) * arma::eye<arma::mat>(3, 3);

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

	this -> original_to_principal_dcm = axes.t();

}


ShapeModel::ShapeModel(std::string ref_frame_name,
                       FrameGraph * frame_graph) {
	this -> frame_graph = frame_graph;
	this -> ref_frame_name = ref_frame_name;
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

std::vector<std::shared_ptr< SBGAT_CORE::Vertex> > * ShapeModel::get_vertices() {
	return &this -> vertices;
}


std::vector<SBGAT_CORE::Facet * > * ShapeModel::get_facets() {
	return &this -> facets;
}

std::vector<std::shared_ptr< SBGAT_CORE::Edge> > * ShapeModel::get_edges() {
	return &this -> edges;
}


void ShapeModel::check_normals_consistency(double tol) const {
	double facet_area_average = 0;

	double sx = 0;
	double sy = 0;
	double sz = 0;

	#pragma omp parallel for reduction(+:facet_area_average,sx,sy,sz) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0; facet_index < this -> facets.size(); ++facet_index) {

		auto facet = this -> facets[facet_index];

		sx += facet -> get_area() * facet -> get_facet_normal() -> at(0);
		sy += facet -> get_area() * facet -> get_facet_normal() -> at(1);
		sz += facet -> get_area() * facet -> get_facet_normal() -> at(2);

		facet_area_average += facet -> get_area();

	}

	arma::vec surface_sum = {sx, sy, sz};

	facet_area_average = facet_area_average / this -> facets.size();
	if (arma::norm(surface_sum) / facet_area_average > tol) {
		throw (std::runtime_error("Normals were incorrectly oriented. norm(sum(n * s))/sum(s):" + std::to_string(arma::norm(surface_sum) / facet_area_average)));
	}

}



void ShapeModel::compute_volume() {
	double volume = 0;

	#pragma omp parallel for reduction(+:volume) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0;
	        facet_index < this -> facets.size();
	        ++facet_index) {

		std::vector<std::shared_ptr<SBGAT_CORE::Vertex > > * vertices = this -> facets[facet_index] -> get_vertices();

		arma::vec * r0 =  vertices -> at(0) -> get_coordinates();
		arma::vec * r1 =  vertices -> at(1) -> get_coordinates();
		arma::vec * r2 =  vertices -> at(2) -> get_coordinates();
		double dv = arma::dot(*r0, arma::cross(*r1 - *r0, *r2 - *r0)) / 6.;
		volume = volume + dv;

	}

	this -> volume = volume;

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


		std::vector<std::shared_ptr<SBGAT_CORE::Vertex > > * vertices = this -> facets[facet_index] -> get_vertices();

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

	arma::vec center_of_mass = {c_x, c_y, c_z};

	this -> cm =  center_of_mass ;

}


void ShapeModel::compute_inertia() {


	double P_xx = 0;
	double P_yy = 0;
	double P_zz = 0;
	double P_xy = 0;
	double P_xz = 0;
	double P_yz = 0;

	double l = std::pow(this -> volume, 1. / 3.);

	#pragma omp parallel for reduction(+:P_xx,P_yy,P_zz,P_xy,P_xz,P_yz) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0;
	        facet_index < this -> facets.size();
	        ++facet_index) {


		std::vector<std::shared_ptr<Vertex > > * vertices = this -> facets[facet_index] -> get_vertices();

		// Normalized coordinates
		arma::vec r0 =  (*vertices -> at(0) -> get_coordinates()) / l;
		arma::vec r1 =  (*vertices -> at(1) -> get_coordinates()) / l;
		arma::vec r2 =  (*vertices -> at(2) -> get_coordinates()) / l;

		double * r0d =  r0. colptr(0);
		double * r1d =  r1. colptr(0);
		double * r2d =  r2. colptr(0);

		double dv = 1. / 6. * arma::dot(r1, arma::cross(r1 - r0, r2 - r0));



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

	// The inertia tensor is centered at the barycenter
	// The parallel axis theorem is used
	// Recall that this inertia tensor is dimensionless


	this -> inertia = I - RBK::tilde(this -> cm / l) * RBK::tilde(this -> cm / l).t();


}


double ShapeModel::get_volume() const {
	return this -> volume;
}


double ShapeModel::get_surface_area() const {
	return this -> surface_area;
}


arma::vec ShapeModel::get_center_of_mass() const {
	return this -> cm;
}

arma::mat ShapeModel::get_original_to_principal_dcm() const {
	return this -> original_to_principal_dcm;
}


void ShapeModel::compute_surface_area() {
	double surface_area = 0;

	#pragma omp parallel for reduction(+:surface_area) if (USE_OMP_SHAPE_MODEL)
	for (unsigned int facet_index = 0; facet_index < this -> facets.size(); ++facet_index) {

		auto facet = this -> facets[facet_index];

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


	std::cout << "xmin : " << xmin << std::endl;
	std::cout << "xmax : " << xmax << std::endl;


	std::cout << "ymin : " << ymin << std::endl;
	std::cout << "ymax : " << ymax << std::endl;


	std::cout << "zmin : " << zmin << std::endl;
	std::cout << "zmax : " << zmax << std::endl;


}


void ShapeModel::set_ref_frame_name(std::string ref_frame_name) {

	this -> ref_frame_name = ref_frame_name;
}

