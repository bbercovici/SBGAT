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

	Edge(std::shared_ptr<Vertex> v0 ,
	     std::shared_ptr<Vertex> v1);

	Vertex * get_v0();
	Vertex * get_v1();
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