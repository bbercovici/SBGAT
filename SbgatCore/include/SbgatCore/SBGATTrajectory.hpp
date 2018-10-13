
#ifndef HEADER_SBGATTRAJECTORY
#define HEADER_SBGATTRAJECTORY

#include <vector>
#include <armadillo>

/**
  @class  SBGATTrajectory
  @author Benjamin Bercovici
  @author Jay McMahon

  @brief  Trajectory-generation class
 
  @details This class can be used to generate trajectories. Trajectories can be loaded from a file
 or generated under a Keplerian dynamics assumption
 
*/
class SBGATTrajectory{

public:

	/**
	Generates a keplerian trajectory at the prescribed time
	@param positions each element in positions holds the inertial cartesian position the corresponding time
	@param velocities each element in velocities holds the inertial cartesian velocity the corresponding time
	@param times vector of times. First time defines the trajectory epoch
	@param elements vector of orbital elements (a,e,i,Omega,omega,M0) with
	- a : semi-major axis (m)
	- e : eccentricity
	- i : inclination in [0,pi] (rad)
	- Omega : right-ascension of ascending node in [0,2pi] (rad)
	- omega : longitude of perigee in [0,2pi] (rad)
	- M0 : mean anomaly at epoch (rad)
	@param mu standard gravitational parameter of orbited body (kg^3/s^2)
	*/
	void GenerateKeplerianTrajectory(
		std::vector<arma::vec> & positions,
		std::vector<arma::vec> & velocities,
		const std::vector<double> &  times,
		const arma::vec & elements,
		const double & mu);

protected:



};




#endif