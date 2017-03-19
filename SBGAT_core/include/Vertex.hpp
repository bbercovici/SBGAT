#ifndef HEADER_VERTEX
#define HEADER_VERTEX
#include <armadillo>
#include "Edge.hpp"
#include "Facet.hpp"

#include <set>

class Edge;
class Facet;

class Vertex {

public:
	arma::vec * get_coordinates() ;
	void set_coordinates(std::shared_ptr<arma::vec> coordinates);

	void add_facet_ownership(Facet * facet);

	std::vector<Facet *>  common_facets(std::shared_ptr<Vertex> vertex) const;
	bool is_owned_by(Facet * facet) const;

	void remove_facet_ownership(Facet * facet);


	unsigned int get_number_of_owning_facets() const ;

protected:
	std::shared_ptr<arma::vec> coordinates;
	std::vector<Facet * > owning_facets;




};


#endif