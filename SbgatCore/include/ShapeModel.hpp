/**
@file   ShapeModel.txt
@Author Benjamin Bercovici (bebe0705@colorado.edu)
@date   May, 2017
@brief  Declaration of the ShapeModel class holding the methods
and pointers representative of a small body shape model
*/



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
Declaration of the ShapeModel class. Represents
a shape parametrized in terms of facets/edges/vertices. In addition, 
edge and facets dyads can be computed for the sake of Polyhedral Gravity Model computation. This class
also provides methods enabling one to compute inertial properties such as the center-of-mass, principal axes and inertia
tensor of the shape being operated on.
*/
class ShapeModel {

public:

	/**
	Constructor
	*/
	ShapeModel();


	/**
	Constructor
	@param frame_graph pointer to the graph storing
	reference frame relationships
	@param frame_graph pointer to the reference frame graph
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
	@param tol numerical tolerance ,i.e value under which the laplacian of the "surface field"
		below which the point is considered outside
	@return true if point is contained inside the shape, false otherwise
	*/
	bool contains(double * point, double tol = 1e-6) ;

	/**
	Checks whether the normals were consistently oriented. If not,
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
	Translates and rotates the shape model so as to have (0,0,0) aligned with its
	barycenter and (1,0,0), (0,1,0) and (0,0,1) with its principal axes
	@param enforce_principal_axes true if the shape must be ligned up 
	with its principal axes of inertia
	*/

	void shift_rotate_to_principal_frame(bool enforce_principal_axes);

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
	Returns the center of mass of the shape expressed in the reference
	frame used before the most recent call to shift_rotate_to_principal_frame
	@return center of mass coordinates in reference frame in use before calling shift_rotate_to_principal_frame
	*/
	arma::vec get_center_of_mass() const;


	/**
	Returns the DCM orienting the shape's principal axes and
	the reference frame used before the most recent call to shift_rotate_to_principal_frame
	@return DCM : P_principal_axis_at_loading_time = DCM * P_in_original_frame
	*/
	arma::mat get_original_to_principal_dcm() const;


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
	Applies a rotation that aligns the body
	with its principal axes.

	This assumes that the body has been shifted so
	that (0,0,0) lies at its barycenter . The resulting inertia tensor is diagonal

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
	Returns the non-dimensional inertia tensor of the body in a
	barycentric frame (rho == 1, l = (volume)^(1/3))
	@return inertia tensor
	*/
	arma::mat get_inertia() const;

protected:

	/**
	Shifts the coordinates of the shape model
	so as to have (0,0,0) aligned with its barycenter
	*/
	void shift_to_barycenter();

	void compute_surface_area();
	void compute_volume();
	void compute_center_of_mass();
	void compute_inertia();
	void compute_principal_axes();

	std::vector<Facet * >  facets;
	std::vector<std::shared_ptr< Edge> >  edges;
	std::vector<std::shared_ptr< Vertex> >  vertices;

	double volume;
	double surface_area;

	bool barycenter_aligned ;
	bool principal_axes_aligned ;


	arma::vec cm;
	arma::mat original_to_principal_dcm;

	arma::mat inertia;


	FrameGraph * frame_graph;
	std::string ref_frame_name;




};
}
#endif