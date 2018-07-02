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

/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SBGATObsLightcurve.cpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <SBGATObsLightcurve.hpp>
#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>
#include <SBGATMassProperties.hpp>
#include <vtkExtractHistogram2D.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <vtkPNGWriter.h>
#include <vtkImageCast.h>
#include <vtkPointData.h>

#include <boost/progress.hpp>



vtkStandardNewMacro(SBGATObsLightcurve);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATObsLightcurve::SBGATObsLightcurve(){

  this -> SetNumberOfOutputPorts(0);
  this -> SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATObsLightcurve::~SBGATObsLightcurve(){

}

int SBGATObsLightcurve::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){

  this -> bspTree_vec.clear();
  this -> polydata_vec.clear();

  // Processing the primary
  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkPolyData * primary = vtkPolyData::SafeDownCast(inInfo0->Get(vtkDataObject::DATA_OBJECT()));
  vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
  tree -> SetDataSet(primary);
  tree -> BuildLocator();
  this -> bspTree_vec.push_back(tree);
  
  this -> polydata_vec.push_back(primary);


  // Processing the secondary, if any
  if(this -> GetNumberOfInputConnections(1) > 0){
   vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);
   vtkPolyData * secondary =  vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));
   vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
   tree -> SetDataSet(secondary);
   tree -> BuildLocator();
   this -> bspTree_vec.push_back(tree);
   this -> polydata_vec.push_back(secondary);
 }

 this -> number_of_bodies = polydata_vec.size();


 return 1;
}



int SBGATObsLightcurve::FillInputPortInformation( int port, vtkInformation* info ){
  if ( port == 0 ){
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;

  }
  else if(port == 1){
   info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
   info->Set(vtkAlgorithm::INPUT_IS_OPTIONAL(), true);
   info->Set(vtkAlgorithm::INPUT_IS_REPEATABLE(), true);
   return 1;

 }

 return 0;
}



void SBGATObsLightcurve::CollectMeasurementsSimpleSpin(
  std::vector<std::array<double, 2> > & measurements,
  const int & N,
  const double & dt,
  const std::vector<double> & period_vec,
  const arma::vec & sun_dir,
  const arma::vec & observer_dir,
  const std::vector<arma::vec> & positions_vec,
  const std::vector<arma::vec> & spin_vec,
  const bool & penalize_indicence){

  if(period_vec.size() != spin_vec.size() || period_vec.size() != positions_vec.size())
    throw(std::runtime_error("Incompatible input dimensions"));


  // Containers
  std::vector<std::vector<int> > facets_in_view;
  std::vector<arma::mat> BN_dcms_vec ;


  // Pre-allocating for all inputs  
  for (int i = 0; i < this -> number_of_bodies; ++i){
    double w = 2 * arma::datum::pi / period_vec[i];
    facets_in_view.push_back(std::vector<int>());
    BN_dcms_vec.push_back(RBK::prv_to_dcm(spin_vec[i] * dt * w ));
  }

  // The surface area of the largest facet 
  double max_area = -1;

  // First, only facets that are in view of the sun
  // and the observer (based on their normal orientation) are kept
  // Then, the facets that are in view of the sun are ray-traced to the sun and to the observer
  // checking with potential interects with all bodies

  std::vector<arma::vec> dir_to_check_vec = {observer_dir,sun_dir};

  for (int i = 0; i < this -> number_of_bodies; ++i){

    this -> prefind_facets_inview(i,
      facets_in_view[i],
      dir_to_check_vec,
      max_area,
      BN_dcms_vec,
      positions_vec);
  }

  for (int i = 0; i < this -> number_of_bodies; ++i){

    std::array<double, 2>  measurements_temp;
    measurements_temp[0] = dt;
    measurements_temp[1] = 0;

    this -> reverse_ray_trace(measurements_temp,
      i,
      facets_in_view,
      sun_dir,
      observer_dir,
      N,
      max_area,
      penalize_indicence,
      BN_dcms_vec,
      positions_vec);

    measurements.push_back(measurements_temp);

  }

}



void SBGATObsLightcurve::CollectMeasurementsArbitrarySpin(
  std::vector<std::array<double, 2> > & measurements,
  const int & N,
  const double & dt,
  const arma::vec & sun_dir,
  const arma::vec & observer_dir,
  const std::vector<arma::vec> & positions_vec,
  const std::vector<arma::mat> & BN_dcms_vec,
  const bool & penalize_indicence){


  if(positions_vec.size() != BN_dcms_vec.size())
    throw(std::runtime_error("Incompatible input dimensions"));
  
  // Containers
  std::vector<std::vector<int> > facets_in_view;

  // Pre-allocating for all inputs  
  for (int i = 0; i < this -> number_of_bodies; ++i){
    facets_in_view.push_back(std::vector<int>());
  }

  // The surface area of the largest facet 
  double max_area = -1;

  // First, only facets that are in view of the sun
  // and the observer (based on their normal orientation) are kept
  // Then, the facets that are in view of the sun are ray-traced to the sun and to the observer
  // checking with potential interects with all bodies
  std::vector<arma::vec> dir_to_check_vec = {observer_dir,sun_dir};
  
  for (int i = 0; i < this -> number_of_bodies; ++i){

    this -> prefind_facets_inview(i,
      facets_in_view[i],
      dir_to_check_vec,
      max_area,
      BN_dcms_vec,
      positions_vec);
  }

  for (int i = 0; i < this -> number_of_bodies; ++i){

    std::array<double, 2>  measurements_temp = {dt, 0};
   
    this -> reverse_ray_trace(measurements_temp,
      i,
      facets_in_view,
      sun_dir,
      observer_dir,
      N,
      max_area,
      penalize_indicence,
      BN_dcms_vec,
      positions_vec);

    measurements.push_back(measurements_temp);

  }

}



void SBGATObsLightcurve::SaveLightCurveData(const std::vector<std::array<double, 2> > & measurements,
  std::string savepath){

  arma::mat time_luminosity_series(measurements.size(),2);

  for (unsigned int i = 0; i < measurements.size(); ++i){
    time_luminosity_series(i, 0) = measurements[i][0];
    time_luminosity_series(i, 1) = measurements[i][1] ;
  }

  time_luminosity_series.save(savepath,arma::raw_ascii);

}



void SBGATObsLightcurve::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATObsLightcurve::PrintTrailer(ostream& os, vtkIndent indent) {

}




//----------------------------------------------------------------------------
void SBGATObsLightcurve::PrintSelf(std::ostream& os, vtkIndent indent){

  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
  if (!input)
  {
    return;
  }

}

void SBGATObsLightcurve::reverse_ray_trace(std::array<double, 2>  & measurements_temp,
  const unsigned int & body_index,
  const std::vector<std::vector<int> > & facets_in_view,
  const arma::vec & sun_dir,
  const arma::vec & observer_dir,
  const int N,
  const double max_area,
  const bool penalize_indicence,
  const std::vector<arma::mat> & BN_dcms_vec,
  const std::vector<arma::vec> & positions_vec){

  vtkPolyData * input = this -> polydata_vec[body_index];

  // Ray-tracing tolerance
  double tol = input -> GetLength()/1E6;

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds -> Allocate(VTK_CELL_SIZE);

  arma::vec target_to_sun_dir_body_frame = BN_dcms_vec[body_index] * sun_dir;
  arma::vec target_to_observer_dir_body_frame = BN_dcms_vec[body_index] * observer_dir;

  // The kept facets are then sampled and reverse ray-traced
  for (unsigned int facet_index = 0; facet_index != facets_in_view[body_index].size(); ++facet_index){

    double p0[3];
    double p1[3];
    double p2[3];

    input -> GetCellPoints(facets_in_view[body_index][facet_index],ptIds);
    input -> GetPoint(ptIds->GetId(0), p0);
    input -> GetPoint(ptIds->GetId(1), p1);
    input -> GetPoint(ptIds->GetId(2), p2);

    arma::vec P0 = {p0[0],p0[1],p0[2]};
    arma::vec P1 = {p1[0],p1[1],p1[2]};
    arma::vec P2 = {p2[0],p2[1],p2[2]};

    arma::vec n = arma::cross(P1 - P0, P2 - P0);

    // The number of points sampled from this facet is determined based on 
    // the relative size of this facet compared to the largest one in all the considered shapes
    int N_samples = int( N * arma::norm(n /2) / max_area);
    
    // only need unit normal vector from here
    n = arma::normalise(arma::cross(P1 - P0, P2 - P0));

    double cosi_sun,cosi_obs ;

    // Computing the ray incidence at impact if needed
    if (penalize_indicence){
      cosi_sun = vtkMath::Dot(n.colptr(0),target_to_sun_dir_body_frame.colptr(0));
      cosi_obs = vtkMath::Dot(n.colptr(0),target_to_observer_dir_body_frame.colptr(0));
    }
    else{
      cosi_sun = 1;
      cosi_obs = 1;
    }

    // A maximum of N points are sampled from this facet
    // #pragma omp parallel for
    for (int i = 0; i < N_samples; ++i){

      bool keep_point = true;

      // A random origin point is uniformly drawn from this facet
      arma::vec random = arma::randu<arma::vec>(2);
      double u = random(0);
      double v = random(1);
      arma::vec origin = (1 - std::sqrt(u)) * P0 + std::sqrt(u) * ( 1 - v ) * P1 + std::sqrt(u) * v * P2 ;

    // Derived points are expressed in the body reference frame
      arma::vec point_above_surface = origin + 3 * tol * n; // the origin of the ray is moved 3*tol above the surface
      arma::vec point_towards_sun = point_above_surface + input -> GetLength() * 1e3 * target_to_sun_dir_body_frame;
      arma::vec point_towards_observer = point_above_surface + input -> GetLength() * 1e3 * target_to_observer_dir_body_frame;

    // Checking if sampled point is in view of sun
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      vtkSmartPointer<vtkPoints> verts = vtkSmartPointer<vtkPoints>::New();

      for (unsigned int considered_body_index = 0; considered_body_index < facets_in_view.size(); ++considered_body_index){

        // The derived points are converted into the considered body's frame
        arma::mat dcm = BN_dcms_vec[considered_body_index] * BN_dcms_vec[body_index].t();
        arma::vec point_above_surface_considered_body = dcm * (point_above_surface  + positions_vec[body_index] - positions_vec[considered_body_index]);
        arma::vec point_towards_sun_considered_body = dcm * (point_towards_sun + positions_vec[body_index] - positions_vec[considered_body_index]);
        arma::vec point_towards_observer_considered_body = dcm * (point_towards_observer + positions_vec[body_index] - positions_vec[considered_body_index]);


        // Checking if the sun view is obscured by the considered body
        this -> bspTree_vec[considered_body_index] -> IntersectWithLine(point_above_surface_considered_body.colptr(0), 
          point_towards_sun_considered_body.colptr(0), tol, verts, cellIds);

        // If there is an intersect, this point is shadowed by the considered body
        if (verts -> GetNumberOfPoints() > 0){
          keep_point = false;
          break;
        }
        

        // Checking if the observer view is obscured
        this -> bspTree_vec[body_index] -> IntersectWithLine(point_above_surface_considered_body.colptr(0), 
          point_towards_observer_considered_body.colptr(0), tol, verts, cellIds);

        // If there is an intersect, this point is shadowed by the considered body
        if (verts -> GetNumberOfPoints() > 0){
          keep_point = false;
          break;
        }
        

      }

      // If this point was not obscured, the return is weighed by the incidence on the inbound and outbout rays

      if (keep_point){
        measurements_temp[1] += cosi_sun * cosi_obs;
      }


    }
  }




}








