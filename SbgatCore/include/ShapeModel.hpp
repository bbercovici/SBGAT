
#ifndef HEADER_SHAPEMODEL
#define HEADER_SHAPEMODEL

#include <string>
#include <iostream>
#include <armadillo>
#include <set>
#include <map>
#include <limits>

#include "Facet.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"
#include "OMP_flags.hpp"
#include "FrameGraph.hpp"


namespace SBGAT_CORE {


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
	Saves the shape model in the form of an .obj file
	@param path Location of the saved file
	*/
	void save(std::string path) const;


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
	Update all the facets of the shape model
	*/
	void update_facets() ;

	/**
	Update all the edges of the shape model
	*/
	void update_edges() ;


	/**
	Updates the specified facets of the shape model. Ensures consistency between the vertex coordinates
	and the facet surface area, normals, dyads and centers.
	@param facets Facets to be updated
	*/
	void update_facets(std::set<Facet *> & facets);


	/**
	Shifts the coordinates of the shape model
	so as to have (0,0,0) aligned with its barycenter

	The resulting barycenter coordinates are (0,0,0)
	*/
	void shift_to_barycenter();

	/**
	Applies a rotation that aligns the body
	with its principal axes.

	This assumes that the body has been shifted so
	that (0,0,0) lies at its barycenter

	The resulting inertia tensor is diagonal

	Undefined behavior if
	the inertia tensor has not been computed beforehand
	*/
	void align_with_principal_axes();

	/**
	Checks if every facets in the shape model has good quality.
	@param min_facet_angle Minimum facet vertex angle indicating degeneracy
	@param min_edge_angle Minimum edge angle indicating degeneracy
	If not, some facets are recycled until the mesh becomes satisfying
	*/
	void enforce_mesh_quality(double min_facet_angle, double min_edge_angle) ;

	/**
	Returns the non-dimensional inertia tensor of the body in the body-fixed
	principal axes. (rho == 1, l = (volume)^(1/3))
	@return principal inertia tensor
	*/
	arma::mat get_inertia() const;




protected:


	void compute_surface_area();
	void compute_volume();
	void compute_center_of_mass();
	void compute_inertia();

	std::vector<Facet * >  facets;
	std::vector<std::shared_ptr< Edge> >  edges;
	std::vector<std::shared_ptr< Vertex> >  vertices;

	double volume;
	double surface_area;

	arma::vec cm;
	arma::mat inertia;


	FrameGraph * frame_graph;
	std::string ref_frame_name;




};
}
#endif