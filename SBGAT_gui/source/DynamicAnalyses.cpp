#include "DynamicAnalyses.hpp"

DynamicAnalyses::DynamicAnalyses(ShapeModel * shape_model) {
	this -> shape_model = shape_model;
}

void DynamicAnalyses::compute_pgm(bool return_pgm ) {

	arma::mat gravity_accelerations = arma::mat(this -> shape_model -> get_NFacets(), 3);

	boost::progress_display progress(this -> shape_model -> get_NFacets()) ;

	for (unsigned int facet = 0; facet < this -> shape_model -> get_NFacets(); ++facet) {

		arma::uvec vertices_in_facet = this -> shape_model -> get_vertex_indices_in_facet(facet);

		arma::vec facet_center = 1. / 3. * (
		                             this -> shape_model -> get_vertex(vertices_in_facet(0) ) +
		                             this -> shape_model -> get_vertex(vertices_in_facet(1) ) +
		                             this -> shape_model -> get_vertex(vertices_in_facet(2) ));

		gravity_accelerations.row(facet) = this -> pgm_acceleration(facet_center).t();
		++progress;
	}

}

arma::vec DynamicAnalyses::pgm_acceleration(arma::vec & point ) const {

	arma::vec acceleration = {0, 0, 0};

	// Facet loop
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

		acceleration += this -> shape_model -> get_F_dyad(facet) * r1 * wf;

	}

	// Edge loop
	for (unsigned int edge_index = 0; edge_index < this -> shape_model -> get_NEdges(); ++ edge_index) {

		std::set<unsigned int> edge = this -> shape_model -> get_edge_from_edge_index(edge_index);

		arma::vec r1 = this -> shape_model -> get_vertex(*edge.begin()) - point;
		arma::vec r2 = this -> shape_model -> get_vertex(*std::next(edge.begin())) - point;

		double R1 = arma::norm(r1);
		double R2 = arma::norm(r2);
		double Re = arma::norm(r2 - r1);
		double Le = std::log((R1 + R2 + Re) / (R1 + R2 - Re));

		acceleration -= this -> shape_model -> get_E_dyad(edge_index) * r1 * Le;

	}


	return acceleration;

}
