#include "DynamicAnalyses.hpp"

DynamicAnalyses::DynamicAnalyses(ShapeModel * shape_model) {
	this -> shape_model = shape_model;
}

void DynamicAnalyses::compute_pgm(double density, bool return_pgm ) {

	arma::mat gravity_accelerations = arma::mat(3, this -> shape_model -> get_NFacets());

	boost::progress_display progress(this -> shape_model -> get_NFacets()) ;

	for (unsigned int facet = 0; facet < this -> shape_model -> get_NFacets(); ++facet) {

		arma::uvec vertices_in_facet = this -> shape_model -> get_vertex_indices_in_facet(facet);

		arma::vec facet_center = 1. / 3. * (
		                             this -> shape_model -> get_vertex(vertices_in_facet(0) ) +
		                             this -> shape_model -> get_vertex(vertices_in_facet(1) ) +
		                             this -> shape_model -> get_vertex(vertices_in_facet(2) ));

		gravity_accelerations.col(facet) = this -> pgm_acceleration(facet_center, density);
		++progress;
	}

}

bool DynamicAnalyses::is_inside(arma::vec point, double tol) const {

	double lagrangian = 0;

	// Facet loop
	#pragma omp parallel for reduction(+:lagrangian)
	for (unsigned int facet = 0; facet < this -> shape_model -> get_NFacets(); ++ facet) {

		arma::uvec vertices_in_facet = this -> shape_model -> get_vertex_indices_in_facet(facet);


		arma::vec r1 = this -> shape_model -> get_vertex(vertices_in_facet(0) ) - point;
		arma::vec r2 = this -> shape_model -> get_vertex(vertices_in_facet(1) ) - point;
		arma::vec r3 = this -> shape_model -> get_vertex(vertices_in_facet(2) ) - point;


		double R1 = arma::norm(r1);
		double R2 = arma::norm(r2);
		double R3 = arma::norm(r3);

		double wf = 2 * std::atan2(arma::dot(r1, arma::cross(r2, r3)),
		                           R1 * R2 * R3 + R1 * arma::dot(r2, r3)
		                           + R2 * arma::dot(r3, r1)
		                           + R3 * arma::dot(r1, r2));

		lagrangian += wf;

	}

	if (std::abs(lagrangian) < tol) {
		return false;
	}
	else {
		return true;
	}

}

arma::vec DynamicAnalyses::pgm_acceleration(arma::vec & point , double density) const {

	double ax = 0;
	double ay = 0;
	double az = 0;

	arma::umat * facet_vertices_pointer =  this -> shape_model -> get_facet_vertices_pointer();
	arma::mat * vertices_pointer = this -> shape_model -> get_vertices_pointer();

	// Facet loop
	#pragma omp parallel for reduction(+:ax,ay,az)
	for (unsigned int facet = 0; facet < this -> shape_model -> get_NFacets(); ++ facet) {

		// arma::uvec vertices_in_facet_indices = this -> shape_model -> get_vertex_indices_in_facet(facet);

		// unsigned int * vertices_in_facet_indices_pointer = this -> shape_model -> get_vertex_indices_in_facet_pointer(facet);



		// double * v0_p = vertices_pointer -> colptr(vertices_in_facet_indices_pointer[0]);
		// double * v1_p = vertices_pointer -> colptr(vertices_in_facet_indices_pointer[1]);
		// double * v2_p = vertices_pointer -> colptr(vertices_in_facet_indices_pointer[2]);

		// arma::vec r1 = vertices_in_facet.col(0) - point;
		// arma::vec r2 = vertices_in_facet.col(1) - point;
		// arma::vec r3 = vertices_in_facet.col(2) - point;

		arma::vec r1 = vertices_pointer -> unsafe_col(
		                   facet_vertices_pointer -> unsafe_col(facet)(0)
		               ) - point;
		arma::vec r2 = vertices_pointer -> unsafe_col(
		                   facet_vertices_pointer -> unsafe_col(facet)(1)
		               ) - point;
		arma::vec r3 = vertices_pointer -> unsafe_col(
		                   facet_vertices_pointer -> unsafe_col(facet)(2)
		               ) - point;

		// arma::vec r1 = *this -> shape_model -> get_vertices_pointer(vertices_in_facet_indices_pointer[0]) - point;
		// arma::vec r2 = *this -> shape_model -> get_vertices_pointer(vertices_in_facet_indices_pointer[1]) - point;
		// arma::vec r3 = *this -> shape_model -> get_vertices_pointer(vertices_in_facet_indices_pointer[2]) - point;

		// arma::vec r1 = vertices_in_facet.col(vertices_in_facet_indices_pointer[0]) - point;
		// arma::vec r2 = vertices_in_facet.col(vertices_in_facet_indices_pointer[1]) - point;
		// arma::vec r3 = vertices_in_facet.col(vertices_in_facet_indices_pointer[2]) - point;

		// arma::vec r1 = vertices_pointer -> col(vertices_in_facet_indices(0)) - point;
		// arma::vec r2 = vertices_pointer -> col(vertices_in_facet_indices(1)) - point;
		// arma::vec r3 = vertices_pointer -> col(vertices_in_facet_indices(2)) - point;


		double R1 = arma::norm(r1);
		double R2 = arma::norm(r2);
		double R3 = arma::norm(r3);

		double wf = 2 * std::atan2(arma::dot(r1, arma::cross(r2, r3)),
		                           R1 * R2 * R3 + R1 * arma::dot(r2, r3)
		                           + R2 * arma::dot(r3, r1)
		                           + R3 * arma::dot(r1, r2));

		arma::vec acc_inc = (this -> shape_model -> get_F_dyad(facet) * r1 * wf);

		ax += acc_inc(0);
		ay += acc_inc(1);
		az += acc_inc(2);
	}



	// Edge loop
	#pragma omp parallel for reduction(-:ax,ay,az)
	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {

		unsigned int v0_index = this -> shape_model -> get_vertex_global_index_from_edge_index(0, edge_index);
		unsigned int v1_index = this -> shape_model -> get_vertex_global_index_from_edge_index(1, edge_index);

		// arma::uvec vertices_in_edge_indices = {v0_index, v1_index};
		// arma::mat vertices_in_edge =  this -> shape_model -> get_vertex_bloc(vertices_in_edge_indices);

		// arma::vec r1 = vertices_in_edge.col(0) - point;
		// arma::vec r2 = vertices_in_edge.col(1) - point;

		arma::vec r1 = vertices_pointer -> unsafe_col(
		                   v0_index) - point;
		arma::vec r2 = vertices_pointer -> unsafe_col(
		                   v1_index) - point;

		double R1 = arma::norm(r1);
		double R2 = arma::norm(r2);
		double Re = arma::norm(r2 - r1);
		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));

		arma::vec acc_dec = this -> shape_model -> get_E_dyad(edge_index) * r1 * Le;

		ax -= acc_dec(0);
		ay -= acc_dec(1);
		az -= acc_dec(2);

	}

	arma::vec acceleration = {ax, ay, az};
	acceleration = acceleration * arma::datum::G * density;

	return acceleration;

}
