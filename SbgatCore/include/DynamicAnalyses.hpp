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