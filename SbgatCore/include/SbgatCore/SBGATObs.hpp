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




#ifndef SBGATOBS_HEADER
#define SBGATOBS_HEADER


#include <vtkModifiedBSPTree.h>
#include <vtkPolydata.h>
#include <armadillo>

class SBGATObs {

public:

	/**
 	 Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters (default)
 	*/
	void SetScaleMeters() { this -> scaleFactor = 1; }

	/**
	Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
	*/
	void SetScaleKiloMeters() { this -> scaleFactor = 1000; }

protected:

	/**
	Determines what facets may be illuminated by the sun. Does not check for intersects and shadowing, 
	only compares the outbound normal orientation to the sun/facet direction. If 
	a facet is found to be visible, its index is added to facets_in_view
	@param body_index index of considered body
	@param facets_in_view reference to vector holding indices of (maybe) illuminated facets
	@param dir_to_check_vec vector of directions (unit vectors) that must point above a facet's horizon for a facet to be considered as visible
	@param max_area maximum area of all facets (illuminated/shaded facets)
	@param BN_dcms_vec vector holding the DCMs orienting the body frame of each body w/r to inertial
	@param positions_vec vector holding the position vector of the CM of each body w/r to the primary
	*/
	void prefind_facets_inview(
		const unsigned int & body_index,
		std::vector<int> & facets_in_view,
		const std::vector<arma::vec> & dir_to_check_vec,
		double & max_area,
		const std::vector<arma::mat> & BN_dcms_vec,
		const std::vector<arma::vec> & positions_vec);


	std::vector<vtkSmartPointer<vtkModifiedBSPTree>> bspTree_vec;

	std::vector<vtkPolyData *> polydata_vec;
	double scaleFactor = 1;
	double max_value;
	int number_of_bodies;

};








#endif