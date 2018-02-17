#ifndef FACETRESULTS_HEADER
#define FACETRESULTS_HEADER

#include <armadillo>


namespace SBGAT_CORE {


/**
Holds pointers to the scalar and vector results held by this facet that
are obtained from analyses
*/
class FacetResults {
public:
	/**
	Constructor
	*/
	FacetResults();

	/**
	Get gravity potential (kg * m^2/s^2)
	@return gravity potential
	*/
	double get_grav_potential() const;

	/**
	Get gravitational slope (deg)
	@return gravity slope
	*/
	double get_grav_slope() const;

	/**
	Get gravity potential (kg * m^2/s^2)
	@param gravity potential
	*/
	void set_grav_potential(double grav_potential);

	/**
	Get gravitational slope (deg)
	@param gravity potential
	*/
	void set_grav_slope(double slope);


	/**
	Return gravity acceleration (kg * m/s^2)
	@return gravity acceleration
	*/
	arma::vec * get_grav_acceleration();

	/**
	Set gravity acceleration (kg * m/s^2)
	@param gravity acceleration
	*/
	void set_grav_acceleration(arma::vec & gravity_acceleration);


protected:
	double grav_potential;
	double grav_slope;

	arma::vec grav_acceleration;





};
}


#endif
