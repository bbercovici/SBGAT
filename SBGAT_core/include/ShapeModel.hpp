
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

class ShapeModel {

public:

	void compute_dyads();
	void compute_F_dyads();
	void compute_F_dyad(unsigned int facet);
	void construct_edges();

	void compute_normals() ;

	void compute_E_dyad(const std::pair < std::set<unsigned int> ,
	                    std::set< unsigned int > > & edge);
	void compute_E_dyads();

	void load(const std::string & filename);
	void save(const std::string & filename) const;

	void load_normals(const std::string & filename);
	void save_normals(const std::string & filename) const ;

	void load_F_dyads(const std::string & filename);
	void save_F_dyads(const std::string & filename) const ;

	void load_E_dyads(const std::string & filename);
	void save_E_dyads(const std::string & filename) const ;

	arma::mat get_F_dyad(unsigned int facet) const ;
	void set_F_dyad(unsigned int facet, arma::mat dyad) ;

	arma::mat get_E_dyad(unsigned int edge_index) const ;
	void set_E_dyad(unsigned int edge_index, arma::mat & dyad) ;

	std::set<unsigned int> get_edge_from_edge_index(unsigned int edge_index) const ;
	unsigned int get_vertex_global_index_from_edge_index(unsigned int vertex_local_index, unsigned int edge_index) const ;

	unsigned int get_NFacets() const ;
	unsigned int get_NVertices() const ;
	unsigned int get_NEdges() const ;


	void set_NFacets(unsigned int nfacet)  ;
	void set_NVertices(unsigned int nvertices)  ;
	void set_NEdges(unsigned int nedges)  ;

	arma::uvec get_vertex_indices_in_facet(unsigned int facet) const;
	arma::vec get_vertex(unsigned int vertex_index) const;
	arma::mat get_vertex_bloc(arma::uvec & cols_to_extract) const ;


	void check_normals_consistency(double tol = 1e-3) const;

	void set_vertices(arma::mat vertices);
	void set_facet_vertices(arma::umat facet_vertices);


protected:

	arma::mat vertices;
	arma::umat facet_vertices;
	arma::mat facet_normals;
	arma::cube E_dyads;
	arma::cube F_dyads;

	arma::umat edges_vertices_indices;

	unsigned int NFacets;
	unsigned int NVertices;
	unsigned int NEdges;

	std::map<std::set<unsigned int> , std::set<unsigned int> > edges_to_facets;
	std::map< std::set<unsigned int> , unsigned int> edges_to_edges_index;





};

#endif