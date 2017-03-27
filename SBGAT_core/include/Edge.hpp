#ifndef HEADER_EDGE
#define HEADER_EDGE
#include <armadillo>
#include "Facet.hpp"
#include "Vertex.hpp"
#include <memory>
#include <set>

class Facet;
class Vertex;

class Edge {

public:

	/**
	Constructor
	@param v0 Pointer to the first edge vertex
	@param v1 Pointer to the second edge vertex
	*/
	Edge(std::shared_ptr<Vertex> v0 ,
	     std::shared_ptr<Vertex> v1);
	/**
	Accessor to the first edge vertex
	@return v0 Pointer to the first edge vertex
	*/
	Vertex * get_v0();

	/**
	Accessor to the second edge vertex
	@return v1 Pointer to the second edge vertex
	*/
	Vertex * get_v1();

	/**
	Returns a pointer to the edge dyad
	@return E Edge dyad
	*/
	arma::mat * get_edge_dyad();


protected:

	void compute_dyad();

	std::shared_ptr<Vertex> v0;
	std::shared_ptr<Vertex> v1;

	Facet * fA;
	Facet * fB;

	std::shared_ptr<arma::mat> edge_dyad;

};
#endif