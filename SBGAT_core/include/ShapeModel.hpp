
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
#include <limits>

#include "Facet.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"
#include "OMP_flags.hpp"
#include "FrameGraph.hpp"


/**
Declaration of the ShapeModel class. Effectively represents
a shape parametrized in terms of facets/edges/vertices. The topology
information is stored in the facets/vertices whereas the edges
store relevant quantities for variational methods such as
the Polyhedron Gravity Model Evaluation
*/
class ShapeModel {

public:

	/**
	Constructor
	*/
	ShapeModel();


	/**
	Constructor
	@param frame_graph Pointer to the graph storing
	reference frame relationships
	@param frame_graph Pointer to the reference frame graph
	*/
	ShapeModel(std::string ref_frame_name,
	           FrameGraph * frame_graph);


	/**
	Destructor
	*/
	~ShapeModel();

	/**
	Returns number of facets
	@return number of facets
	*/
	unsigned int get_NFacets() const ;

	/**
	Returns number of vertices
	@return number of vertices
	*/
	unsigned int get_NVertices() const ;

	/**
	Returns number of edges
	@return number of edges
	*/
	unsigned int get_NEdges() const ;


	/**
	Determines whether the provided point lies inside or outside the shape model.
	The shape model must have a closed surface for this method to be trusted
	@param point coordinates of the point to be tested expressed in the shape model frame
	@param tol numerical tolerance ,i.e value under which the lagrangian of the "surface field"
		below which the point is considered outside
	@return true if point is contained inside the shape, false otherwise
	*/
	bool contains(double * point, double tol = 1e-6) ;


	/**
	Checks that the normals were consistently oriented. If not,
	the ordering of the vertices in the provided shape model file is incorrect
	@param tol numerical tolerance (if consistent: norm(Sum(oriented_surface_area)) / average_facet_surface_area << tol)
	*/
	void check_normals_consistency(double tol = 1e-3) const;


	/**
	Augment the internal container storing facets with a new (and not already inserted)
	one
	@param facet pointer to the new facet to be inserted
	*/
	void add_facet(Facet * facet);

	/**
	Augment the internal container storing edge with a new (and not already inserted)
	one
	@param edge pointer to the new edge to be inserted
	*/
	void add_edge(std::shared_ptr<Edge> edge);

	/**
	Augment the internal container storing vertices with a new (and not already inserted)
	one
	@param vertex pointer to the new vertex to be inserted
	*/
	void add_vertex(std::shared_ptr<Vertex> vertex);

	/**
	Defines the reference frame attached to the shape model
	@param ref_frame Pointer to the reference frame attached
	to the shape model
	*/
	void set_ref_frame_name(std::string ref_frame_name);

	/**
	Pointer to the shape model's vertices
	@return vertices pointer to the vertices
	*/
	std::vector<std::shared_ptr< Vertex> > * get_vertices();

	/**
	Pointer to the shape model's facets
	@return facets pointer to the facets
	*/
	std::vector<Facet * >  * get_facets();

	/**
	Pointer to the shape model's edges
	@return edges pointer to the edges
	*/
	std::vector<std::shared_ptr< Edge> > * get_edges();

	/**
	Returns the dimensions of the bounding box
	@param Bounding box dimension to be computed (xmin,ymin,zmin,xmax,ymax,zmax)
	*/
	void get_bounding_box(double * bounding_box) const;

	/**
	Returns the surface area of the shape model
	@return surface area (U^2 where U is the unit of the shape coordinates)
	*/
	double get_surface_area() const;

	/**
	Returns the volume of the provided shape model
	@return volume (U^2 where U is the unit of the shape coordinates)
	*/
	double get_volume() const;

	/**
	Returns the location of the center of mass
	@return pointer to center of mass
	*/
	arma::vec * get_center_of_mass();

	/**
	Returns the name of the reference frame attached to this
	ref frame
	@return name of reference frame
	*/
	std::string get_ref_frame_name() const;


	/**
	Updates the values of the center of mass, volume, surface area
	*/
	void update_mass_properties();

	/**
	Shifts the coordinates of the shape model
	@param x Shifting to apply to the coordinates
	*/
	void shift(arma::vec x);


protected:


	void compute_surface_area();
	void compute_volume();
	void compute_center_of_mass();

	std::vector<Facet * >  facets;
	std::vector<std::shared_ptr< Edge> >  edges;
	std::vector<std::shared_ptr< Vertex> >  vertices;

	double volume;
	arma::vec cm;
	double surface_area;

	FrameGraph * frame_graph;
	std::string ref_frame_name;


};

#endif