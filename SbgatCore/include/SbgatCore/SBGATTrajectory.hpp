/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/



#ifndef HEADER_SBGATTRAJECTORY
#define HEADER_SBGATTRAJECTORY

#include <vector>
#include <armadillo>

/**
 * @class  SBGATTrajectory
 * @author Benjamin Bercovici
 * @brief  Trajectory-generation class
 *
 * @details This class can be used to generate trajectories. Trajectories can be loaded from a file
 or generated under a Keplerian dynamics assumption
 *
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