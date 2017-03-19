
#ifndef HEADER_SHAPEMODEL
#define HEADER_SHAPEMODEL

#include <string>
#include <iostream>
#include <armadillo>
// #include <assimp/Importer.hpp>
// #include <assimp/scene.h>
// #include <assimp/postprocess.h>
#include <set>
#include <map>

#include "Facet.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"

class ShapeModel {

public:

	~ShapeModel();

	void compute_dyads();
	void compute_F_dyads();
	void compute_F_dyad(unsigned int facet);
	void construct_edges();

	void compute_normals() ;

	void compute_E_dyad(const std::pair < std::set<unsigned int> ,
	                    std::set< unsigned int > > & edge);
	void compute_E_dyads();


	std::shared_ptr<arma::mat> get_F_dyad_ptr(unsigned int facet) ;
	std::shared_ptr<arma::mat> get_E_dyad_ptr(unsigned int edge_index) ;

	std::set<unsigned int> get_edge_from_edge_index(unsigned int edge_index) const ;
	unsigned long long * get_vertex_global_index_from_edge_index(unsigned int edge_index) ;

	unsigned int get_NFacets() const ;
	unsigned int get_NVertices() const ;
	unsigned int get_NEdges() const ;

	bool contains(double * point, double tol = 1e-6) ;

	double * get_facet_normal(unsigned int facet);


	void set_NFacets(unsigned int nfacet)  ;
	void set_NVertices(unsigned int nvertices)  ;
	void set_NEdges(unsigned int nedges)  ;

	arma::uvec get_vertex_indices_in_facet(unsigned int facet) const;
	arma::vec get_vertex(unsigned int vertex_index) const;
	arma::mat get_vertex_bloc(arma::uvec & cols_to_extract) const ;

	unsigned int * get_vertex_indices_in_facet_pointer(unsigned int facet) ;
	arma::mat * get_vertices_pointer() ;
	arma::umat * get_facet_vertices_pointer();

	void check_normals_consistency(double tol = 1e-3) const;


	void add_facet(Facet * facet);
	void add_edge(std::shared_ptr<Edge> edge);
	void add_vertex(std::shared_ptr<Vertex> vertex);

	std::vector<std::shared_ptr< Vertex> > * get_vertices();
	std::vector<Facet * >  * get_facets();
	std::vector<std::shared_ptr< Edge> > * get_edges();


protected:

	std::vector<Facet * >  facets;
	std::vector<std::shared_ptr< Edge> >  edges;
	std::vector<std::shared_ptr< Vertex> >  vertices;

	std::vector< std::shared_ptr< std::vector<Facet * > > > vertices_to_facets;




};

#endif