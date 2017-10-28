/**
@file   DynamicAnalyses.txt
@Author Benjamin Bercovici (bebe0705@colorado.edu)
@date   May, 2017
@brief  Declaration of the DynamicAnalyses class holding the method and members required to perform
a variety of analyses over small body shape models.
*/





#ifndef HEADER_DYNAMICANALYSES
#define HEADER_DYNAMICANALYSES

#include "ShapeModel.hpp"
#include <armadillo>
#include <boost/progress.hpp>
#include "omp.h"
#include "OMP_flags.hpp"


namespace SBGAT_CORE {


	/**
	Declaration of the DynamicAnalyses class. Holds
	methods and members enabling geophysical analyses by
	operating over small body shape models. 
	*/
	class DynamicAnalyses {

	public:
	/**
	Constructor. Creates an instance of a DynamicalAnalyses object
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
	Compute the exterior gravity spherical harmonics coefficients assuming a constant density polyhedron. 
	@param Cnm_total computed Cnm coefficients
	@param Snm_total computed Snm coefficients
	@param n_degree maximum degree of the spherical harmonics expansion
	@param ref_radius normalizing radius (L)
	@param density density of polyhedron (M/L^3)
	@param normalized if true, computed normalized coefficients
	*/
		void compute_exterior_sh_coefs(
			arma::mat & Cnm_total,
			arma::mat & Snm_total,
			unsigned int n_degree,
			double ref_radius,
			double density,
			bool normalized);

	/**
	Computes gravitational slopes, evaluated and stored at the center of each facet.
	compute_pgm_accelerations() must have been called prior to compute_gravity_slopes().
	@param spin_axis Spin axis in body-fixed frame
	@param spin_rate Spin rate (rad/s)
	@return Vector storing the [min,mean,max] slopes (deg)
	*/
		arma::vec compute_gravity_slopes(
			arma::vec spin_axis,
			double spin_rate);

	/**
	Evaluates the acceleration due to gravity at the provided point using a Polyhedron Gravity Model
	@param point body-frame coordinates at which the acceleration is evaluated (L)
	@param density density of the shape model (M/L^3)
	@return PGM acceleration
	*/
		arma::vec pgm_acceleration(double * point , double density) const ;

	/**
	Evaluates the gravity potential at the provided point using a Polyhedron Gravity Model
	@param point body-frame coordinates at which the acceleration is evaluated (L)
	@param density Density of the shape model (M/L^3)
	@return PGM potential
	*/
		double pgm_potential(double * point , double density) const;


	/**
	Evaluates the gravity acceleration from the provided spherical harmonics 
	coefficients. This function should not be used if 
	the query point lies within the circumscribing sphere surrounding the body .

	The query point coordinates must be expressed in the same frame as the one that 
	was used to evaluate the spherical harmonics coefficients 

	This method was reimplemented by Benjamin Bercovici from the original works
	of Yu Takahashi and Siamak Hesar from CU Boulder.

	For more information, one can read the following references

	- 1: S. V. Bettadpur, "Hotine's geopotential formulation: revisited", Bulletin Geodesique (1995) 69:i35-142
  	- 2: R. A. Werner, "Evaluating Descent and Ascent Trajectories Near Non-Spherical Bodies", Technical Support Package
	- 3: L. E. Cunningham, "On the computation of the spherical harmonic terms needed during the numerical integration of the orbital motion of an artificial satellite"

	@param n_degree degree of the expansion
	@param ref_radius reference radius used in the expansion [L]
	@param mu standard gravitational parameter of the attracting body [L^3/s^2]
	@param pos query point position expressed in the body-fixed frame of the 
	attracting body [L]
	@param Cbar matrix of normalized C coefficients
	@param Sbar matrix of normalized S coefficients
	@return spherical harmonics acceleration evaluated in the body-fixed frame
	*/
		arma::vec spherical_harmo_acc(const unsigned int n_degree,
			const double ref_radius,
			const double  mu,
			const arma::vec & pos, 
			const arma::mat & Cbar,
			const arma::mat & Sbar);

	protected:

		ShapeModel * shape_model;

		void GetBnmNormalizedExterior(int n_degree,
			arma::mat & b_bar_real,
			arma::mat & b_bar_imag,
			const arma::vec & pos,
			double ref_radius);




	};

}


#endif