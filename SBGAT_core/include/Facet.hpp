#ifndef HEADER_FACET
#define HEADER_FACET
#include <armadillo>
#include "Edge.hpp"
#include "Vertex.hpp"
#include <memory>

#include <set>

class Edge;
class Vertex;

class Facet {

public:

	/**
	Constructor
	@param vertices pointer to vector storing the vertices owned by this facet
	*/
	Facet(std::shared_ptr< std::vector<std::shared_ptr<Vertex > > > vertices);

	/**
	Not implemented
	*/
	void add_neighbor(std::shared_ptr< Facet > facet);


	/**
	Not implemented
	*/
	void add_edge(std::shared_ptr< Edge > edge);


	/**
	Not implemented
	*/
	void remove_neighbor(std::shared_ptr< Facet > facet);


	/**
	Not implemented
	*/
	void remove_edge(std::shared_ptr< Facet > edge);



	/**
	Get outbound facet normal
	@return pointer to facet normal
	*/
	arma::vec * get_facet_normal()  ;

	/**
	Get facet dyad
	@return pointer to facet dyad
	*/
	arma::mat * get_facet_dyad() ;

	/**
	Not implemented
	*/
	std::set<std::shared_ptr< Facet > > get_neighbors() const;

	/**
	Not implemented
	*/
	std::vector<std::shared_ptr <Edge> > get_facet_edges() const;


	/**
	Returns pointer to the first vertex owned by $this that is 
	neither $v0 and $v1. When $v0 and $v1 are on the same edge, 
	this method returns a pointer to the vertex of $this that is not
	on the edge but still owned by $this
	@param v0 Pointer to first vertex to exclude
	@param v1 Pointer to first vertex to exclude
	@return Pointer to the first vertex of $this that is neither $v0 and $v1
	*/
	std::shared_ptr<Vertex> vertex_not_on_edge(std::shared_ptr<Vertex> v0,
	        std::shared_ptr<Vertex>v1) const ;

	/**
	Not implemented
	*/
	void update();

	/**
	Return facet surface area
	@return facet surface area
	*/
	double get_area() const;

	/**
	Return pointer to facet center
	@return pointer to facet center
	*/
	arma::vec * get_facet_center() ;

	/**
	Return the vector storing the vertices owned by this facet
	@return vector storing the vertices owned by this facet
	*/
	std::vector<std::shared_ptr<Vertex > > * get_vertices() ;



protected:

	void compute_facet_dyad();
	void compute_normal();
	void compute_area();
	void compute_facet_center();


	std::shared_ptr< std::vector<std::shared_ptr<Vertex > > > vertices ;
	std::shared_ptr<arma::mat> facet_dyad;
	std::shared_ptr<arma::vec> facet_normal;
	std::shared_ptr<arma::vec> facet_center;

	std::set<std::shared_ptr< Facet > > neighbors;
	std::vector<std::shared_ptr <Edge> > facet_edges;
	double area;

};
#endif