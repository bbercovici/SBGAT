#ifndef HEADER_DYNAMICANALYSES
#define HEADER_DYNAMICANALYSES

#include "ShapeModel.hpp"
#include <armadillo>
#include <boost/progress.hpp>
#include "omp.h"
#include "OMP_flags.hpp"


namespace SBGAT_CORE {

class DynamicAnalyses {

public:
	/**
	Constructor
	Creates an instance of a DynamicalAnalyses object
	@param shape_model Pointer to shape model for which analysis must be conducted
	*/
	DynamicAnalyses(ShapeModel * shape_model);

	/**
	Evaluates the Polyhedron Gravity Model acceleration (PGM) at the center of each facet of the shape model
	@param density Density of the shape model (kg/m^3)
	*/
	void compute_pgm_accelerations(double density);

	/**
	Evaluates the Polyhedron Gravity Model potentials (PGM) at the center of each facet of the shape model
	@param density Density of the shape model (kg/m^3)
	*/
	void compute_pgm_potentials(double density);


	/**
	Compute the normalized exterior gravity spherical harmonics coefficients assuming a constant density polyhedron
	@param Cnm_total computed Cnm and Snm coefficients
	@param Snm_total computed Snm coefficients
	@param n_degree maximum degree of the spherical harmonics expansion
	@param ref_radius normalizing radius (L)
	@param density density of polyhedron (M/L^3)
	@param normalized if true, returns normalized coefficients
	*/
	void compute_exterior_sh_coefs_normalized(
	    arma::mat & Cnm_total,
	    arma::mat & Snm_total,
	    unsigned int n_degree,
	    double ref_radius,
	    double density,
	    bool normalized);

	/**
	Computes gravitational slopes, evaluated and stored at the center of each facet
	@param spin_axis Spin axis in body-fixed frame
	@param spin_rate Spin rate (rad/s)
	@return Vector storing the [min,mean,max] slopes (deg)
	*/
	arma::vec compute_gravity_slopes(
	    arma::vec spin_axis,
	    double spin_rate);

	/**
	Evaluates the acceleration due to gravity at the provided point using a Polyhedron Gravity Model
	@param point Array of body-frame coordinates at which the acceleration is evaluated
	@param density Density of the shape model (kg/m^3)
	@return PGM acceleration
	*/
	arma::vec pgm_acceleration(double * point , double density) const ;

	/**
	Evaluates the gravity potential at the provided point using a Polyhedron Gravity Model
	@param point Array of body-frame coordinates at which the acceleration is evaluated
	@param density Density of the shape model (kg/m^3)
	@return PGM potential
	*/
	double pgm_potential(double * point , double density) const;

protected:

	ShapeModel * shape_model;

};

}


#endif