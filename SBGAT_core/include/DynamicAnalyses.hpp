#ifndef HEADER_DYNAMICANALYSES
#define HEADER_DYNAMICANALYSES

#include "ShapeModel.hpp"
#include <armadillo>
#include <boost/progress.hpp>
#include "omp.h"
#include "OMP_flags.hpp"

class DynamicAnalyses {

public:
	/**
	Constructor
	Creates an instance of a DynamicalAnalyses object
	@param shape_model Pointer to shape model for which analysis must be conducted
	*/
	DynamicAnalyses(ShapeModel * shape_model);

	/**
	Evaluates the Polyhedron Gravity Model (PGM) at the center of each facet of the shape model
	@param density Density of the shape model (kg/m^3)
	@param XXXXXX
	*/
	void compute_pgm(double density, bool return_pgm = false);

	/**
	Saves the gravity acceleration to file.
	@param path Save path
	*/
	void save_gravity_accelerations_to_path(std::string path) const;

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
	Evaluates the gravity acceleration due to gravity at the provided point using a Polyhedron Gravity Model
	@param point Array of body-frame coordinates at which the acceleration is evaluated
	@param density Density of the shape model (kg/m^3)
	@return PGM acceleration
	*/
	arma::vec pgm_acceleration(double * point , double density) const ;

protected:


	ShapeModel * shape_model;
	arma::mat gravity_accelerations;
	arma::vec slopes;



};


#endif