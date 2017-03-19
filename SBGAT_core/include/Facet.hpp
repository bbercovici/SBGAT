#ifndef HEADER_FACET
#define HEADER_FACET
#include <armadillo>
#include "Edge.hpp"
#include "Vertex.hpp"

#include <set>

class Edge;
class Vertex;

class Facet {

public:

	Facet(std::shared_ptr< std::vector<std::shared_ptr<Vertex > > > vertices);

	void add_neighbor(std::shared_ptr< Facet > facet);
	void add_edge(std::shared_ptr< Edge > edge);

	void remove_neighbor(std::shared_ptr< Facet > facet);
	void remove_edge(std::shared_ptr< Facet > edge);

	arma::vec * get_facet_normal()  ;
	arma::mat * get_facet_dyad() ;

	std::set<std::shared_ptr< Facet > > get_neighbors() const;
	std::vector<std::shared_ptr <Edge> > get_facet_edges() const;

	std::shared_ptr<Vertex> vertex_not_on_edge(std::shared_ptr<Vertex> v0,
	        std::shared_ptr<Vertex>v1) const ;


	void update();
	double get_area() const;
	arma::vec * get_facet_center() ;

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