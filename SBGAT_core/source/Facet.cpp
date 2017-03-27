#include "../include/Facet.hpp"
#include <memory>

Facet::Facet(std::shared_ptr< std::vector<std::shared_ptr<Vertex > > >   vertices) {
	this -> vertices = vertices;

	for (unsigned int vertex_index = 0; vertex_index < this -> vertices -> size(); ++vertex_index) {
		this -> vertices -> at(vertex_index) -> add_facet_ownership(this);
	}

	// Allocating memory for the facet normal
	this -> facet_normal = std::make_shared<arma::vec>(arma::zeros(3));
	this -> compute_normal();

	// Allocating memory for the facet dyad
	this -> facet_dyad = std::make_shared<arma::mat>(arma::zeros(3, 3));
	this -> compute_facet_dyad();

	// Allocating memory for the facet center
	this -> facet_center = std::make_shared<arma::vec>(arma::zeros(3));;
	this -> compute_facet_center();


	// Computing surface area
	this -> compute_area();

}

void Facet::compute_normal() {

	arma::vec * P0 = this -> vertices -> at(0) -> get_coordinates();
	arma::vec * P1 = this -> vertices -> at(1) -> get_coordinates();
	arma::vec * P2 = this -> vertices -> at(2) -> get_coordinates();

	*this -> facet_normal = arma::cross(*P1 - *P0, *P2 - *P0) / arma::norm(arma::cross(*P1 - *P0, *P2 - *P0));
}

void Facet::compute_facet_dyad() {

	*this -> facet_dyad = *this -> facet_normal * (*this -> facet_normal).t();
}

arma::vec * Facet::get_facet_normal()  {
	return this -> facet_normal.get();
}


arma::mat * Facet::get_facet_dyad()  {
	return this -> facet_dyad.get();
}

std::shared_ptr<Vertex> Facet::vertex_not_on_edge(std::shared_ptr<Vertex> v0,
        std::shared_ptr<Vertex>v1) const {
	for (unsigned int i = 0; i < this -> vertices -> size(); ++i) {

		if (this -> vertices -> at(i) != v0 && this -> vertices -> at(i) != v1 ) {
			return this -> vertices -> at(i);
		}
	}
	return nullptr;

}


arma::vec * Facet::get_facet_center()  {
	return this -> facet_center.get();

}

void Facet::compute_facet_center() {

	arma::vec facet_center = arma::zeros(3);

	for (unsigned int vertex_index = 0; vertex_index < this -> vertices -> size(); ++vertex_index) {

		facet_center += *this -> vertices -> at(vertex_index) -> get_coordinates();

	}

	*this -> facet_center = facet_center / this -> vertices -> size();

}

std::vector<std::shared_ptr<Vertex > >  * Facet::get_vertices() {
	return this -> vertices.get();
}


void Facet::compute_area() {
	arma::vec * P0 = this -> vertices -> at(0) -> get_coordinates() ;
	arma::vec * P1 = this -> vertices -> at(1) -> get_coordinates() ;
	arma::vec * P2 = this -> vertices -> at(2) -> get_coordinates() ;
	this -> area = arma::norm( arma::cross(*P1 - *P0, *P2 - *P0)) / 2;
}

double Facet::get_area() const {
	return this -> area;
}