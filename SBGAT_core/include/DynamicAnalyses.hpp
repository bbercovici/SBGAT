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
	Computes gravitational slopes
	@param spin_axis Spin axis in body-fixed frame
	@param spin_rate Spin rate (rad/s)
	*/
	void compute_slopes(
	    arma::vec spin_axis,
	    double spin_rate);

	/**
	Saves the gravitional slopes to file
	@param path Save path
	*/
	void save_slopes(std::string path) const;


	/**
	Saves the gravitional potentials to file
	@param path Save path
	*/
	void save_pgm_potentials(std::string path) const;

	/**
	Saves the gravity acceleration to file.
	@param path Save path
	*/
	void save_pgm_accelerations(std::string path) const;


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
	arma::mat pgm_accelerations;
	arma::vec pgm_potentials;
	arma::vec slopes;



};

}


#endif