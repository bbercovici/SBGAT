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
#include <vtkPolyDataAlgorithm.h>
#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <armadillo>



class VTKFILTERSCORE_EXPORT SBGATObs : public vtkPolyDataAlgorithm {

public:


	static SBGATObs * New();

	vtkTypeMacro(SBGATObs,vtkPolyDataAlgorithm);

	void PrintSelf(std::ostream& os, vtkIndent indent) override;
	void PrintHeader(std::ostream& os, vtkIndent indent) override;
	void PrintTrailer(std::ostream& os, vtkIndent indent) override;

	/**
 	Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters (default)
 	*/
	void SetScaleMeters() { this -> scaleFactor = 1; }

	/**
	Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
	*/
	void SetScaleKiloMeters() { this -> scaleFactor = 1000; }

	/**
	Returns the number of considered bodies
	@return number of considered bodies
	*/
	int get_number_of_bodies() const {return this -> number_of_bodies;}

protected:

	SBGATObs();
	~SBGATObs() override;


	int RequestData(vtkInformation* request,
		vtkInformationVector** inputVector,
		vtkInformationVector* outputVector) override;


	int FillInputPortInformation( int port, vtkInformation* info ) VTK_OVERRIDE;


	/*
	Determines what facets of a given shape satisfy a set of visibility conditions, formulated in terms of a 
	direction that must point above the facet's horizon. Does not check for intersects and shadowing, 
	only compares the outbound normal orientation to the unit vector direction of each condition. If 
	a facet is found to be visible (i.e all conditions are satisfied), 
	its index is added to facets_in_view.
	@param facets_in_view reference to vector holding indices of (maybe) illuminated facets
	@param body_index index of considered body
	@param dir_to_check_vec vector of directions (unit vectors) that must point above a facet's horizon for a facet to be considered as visible
	@param BN_dcms_vec vector holding the DCMs orienting the body frame of each body w/r to inertial
	@param positions_vec vector holding the position vector of the CM of each body w/r to the primary
	*/
	void prefind_facets_inview(
		std::vector<int> & facets_in_view,
		const unsigned int & body_index,
		const std::vector<arma::vec> & dir_to_check_vec,
		const std::vector<arma::mat> & BN_dcms_vec,
		const std::vector<arma::vec> & positions_vec);

	/**
	Computes the largest facet surface area of all of the considered shapes
	@return surface area of largest facet in all of the considered shapes
	*/
	void find_min_facet_surface_area();

	/**
	Checks if the line spanned between provided points intersects
	with any of the considered bodies. 
	@param origin_body_index index of the body from which originates the line to be tested for intersect
	@param start_point_origin_body coordinates of the origin of the line, expressed in the original body frame
	@param end_point_inertial coordinates of the end of the line, expressed in the inertial frame
	@param BN_dcms_vec vector holding the DCMs orienting the body frame of each body w/r to inertial
	@param positions_vec vector holding the position vector of the CM of each body w/r to the primary
	@return true if the line intersects with any of the considered bodies, false otherwise
	*/
	bool check_line_for_intersect(const int & origin_body_index,
		const arma::vec & start_point_origin_body,
		const arma::vec & end_point_inertial,
		const std::vector<arma::mat> & BN_dcms_vec,
		const std::vector<arma::vec> & positions_vec,
		const double & tol) const;


	std::vector<vtkSmartPointer<vtkModifiedBSPTree>> bspTree_vec;
	std::vector<vtkPolyData *> polydata_vec;
	std::vector<arma::vec> center_of_mass_vec;
	
	double scaleFactor = 1;
	double min_area;
	int number_of_bodies;

private:
	SBGATObs(const SBGATObs&) = delete;
	void operator=(const SBGATObs&) = delete;

};








#endif