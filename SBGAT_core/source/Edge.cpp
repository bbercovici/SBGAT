#include "../include/Edge.hpp"

Edge::Edge(std::shared_ptr<Vertex> v0 ,
           std::shared_ptr<Vertex> v1) {

	this ->  v0 = v0;
	this ->  v1 = v1;

	std::vector<Facet *>  common_facets = v0 -> common_facets(v1);

	this -> fA = common_facets.at(0);
	this -> fB = common_facets.at(1);

	this -> edge_dyad = std::make_shared<arma::mat>(arma::zeros(3, 3));


	this -> compute_dyad();

}

Vertex * Edge::get_v0() {
	return this -> v0.get();
}

Vertex * Edge::get_v1() {
	return this -> v1.get();
}

arma::mat * Edge::get_edge_dyad()  {

	return this -> edge_dyad.get();
}


void Edge::compute_dyad() {

	arma::vec * facet_normal_A = this -> fA -> get_facet_normal();
	arma::vec * facet_normal_B =  this -> fB -> get_facet_normal();

	arma::vec edge_dir = arma::normalise(arma::cross(*facet_normal_A, *facet_normal_B) );
	arma::vec edge_normal_A_to_B = arma::normalise(arma::cross(edge_dir, *facet_normal_A));
	arma::vec edge_normal_B_to_A = arma::normalise(arma::cross(edge_dir, *facet_normal_B));

	arma::vec * vertex_in_A_not_on_edge = this -> fA -> vertex_not_on_edge(this -> v0,
	                                      this -> v1) -> get_coordinates();
	arma::vec * vertex_in_B_not_on_edge = this -> fB -> vertex_not_on_edge(this -> v0,
	                                      this -> v1) -> get_coordinates();


	// Consistency check
	if (arma::dot(edge_normal_A_to_B,
	              *this -> v0 -> get_coordinates() - *vertex_in_A_not_on_edge) < 0) {
		edge_normal_A_to_B = - edge_normal_A_to_B;
	}


	if (arma::dot(edge_normal_B_to_A,
	              *this -> v0 -> get_coordinates() - *vertex_in_B_not_on_edge) < 0) {
		edge_normal_B_to_A = - edge_normal_B_to_A;
	}



	*this -> edge_dyad = (*facet_normal_A) * edge_normal_A_to_B.t() + (*facet_normal_B) * edge_normal_B_to_A.t();;

}