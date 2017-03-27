#ifndef HEADER_VERTEX
#define HEADER_VERTEX
#include <armadillo>
#include "Edge.hpp"
#include "Facet.hpp"
#include <memory>

#include <set>

class Edge;
class Facet;

class Vertex {

public:

	/**
	Getter to the vertex's coordinates
	@return coordinates vertex coordinates
	*/
	arma::vec * get_coordinates() ;

	/**
	Setter to the vertex's coordinates
	@param coordinates vertex coordinates
	*/
	void set_coordinates(std::shared_ptr<arma::vec> coordinates);


	/**
	Adds $facet to the vector of Facet * that own this vertex.
	Nothing happens if the facet is already listed
	@param facet Pointer to the facet owning this vertex
	*/
	void add_facet_ownership(Facet * facet);



	/**
	Finds the facets owming both $this and $vertex
	@return commons_facets Vector of Facet * owning the two vertices
	*/
	std::vector<Facet *>  common_facets(std::shared_ptr<Vertex> vertex) const;

	/**
	Determines if $this is owned by $facet
	@param facet Facet whose relationship with the facet is to be tested
	@return true is $this is owned by $facet, false otherwise
	*/
	bool is_owned_by(Facet * facet) const;


	/**
	Delete $facet from the list of Facet * owning $this
	Nothing happens if the facet was not listed (maybe throw a warning)>
	@param facet Pointer to the facet owning this vertex
	*/
	void remove_facet_ownership(Facet * facet);


	/**
	Returns the number of facets owning this vertex
	@return N number of owning of facets
	*/
	unsigned int get_number_of_owning_facets() const ;

protected:
	std::shared_ptr<arma::vec> coordinates;
	std::vector<Facet * > owning_facets;




};


#endif