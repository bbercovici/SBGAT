/*=========================================================================

  Program:   Visualization Toolkit
  Module:    SBGATObsLightcurve.cpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <SBGATObsLightcurve.hpp>
#include <SBGATObs.hpp>

#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
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


SBGATObsLightcurve::SBGATObsLightcurve() : SBGATObs(){

}


SBGATObsLightcurve::~SBGATObsLightcurve(){

}


void SBGATObsLightcurve::CollectMeasurements(
  std::vector<std::array<double, 2> > & measurements,
    const double & time,
    const int & N,
    const arma::vec & sun_dir,
    const arma::vec & observer_dir,
    const std::vector<arma::vec> & positions_vec,
    const std::vector<arma::vec> & velocities_vec,
    const std::vector<arma::vec> & mrps_vec,
    const std::vector<arma::vec> & omegas_vec,
    const bool & penalize_indicence){

 
  if(positions_vec.size() != omegas_vec.size() || positions_vec.size() != velocities_vec.size()|| positions_vec.size() != mrps_vec.size())
    throw(std::runtime_error("Incompatible input dimensions"));


// Containers
  std::vector<std::vector<int> > facets_in_view;
  std::vector<arma::mat> BN_dcms_vec ;



  // Pre-allocating for all inputs  
  for (int i = 0; i < this -> number_of_bodies; ++i){
    facets_in_view.push_back(std::vector<int>());
    BN_dcms_vec.push_back(RBK::mrp_to_dcm(mrps_vec[i]));
  }



  // First, only facets that are in view of the sun
  // and the observer (based on their normal orientation) are kept
  // Then, the facets that are in view of the sun are ray-traced to the sun and to the observer
  // checking with potential interects with all bodies

  std::vector<arma::vec> dir_to_check_vec = {observer_dir,sun_dir};


  for (int i = 0; i < this -> number_of_bodies; ++i){

    this -> prefind_facets_inview(facets_in_view[i],
      i,
      dir_to_check_vec,
      BN_dcms_vec,
      positions_vec);
  }


  std::array<double, 2>  measurements_temp;
  measurements_temp[0] = time;
  measurements_temp[1] = 0;

  this -> reverse_ray_trace(measurements_temp,
    facets_in_view,
    sun_dir,
    observer_dir,
    N,
    penalize_indicence,
    BN_dcms_vec,
    positions_vec);


  measurements.push_back(measurements_temp);


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
  const std::vector<std::vector<int> > & facets_in_view,
  const arma::vec & sun_dir,
  const arma::vec & observer_dir,
  const int N,
  const bool penalize_indicence,
  const std::vector<arma::mat> & BN_dcms_vec,
  const std::vector<arma::vec> & positions_vec){


  

  // The sun is positionned with respect to the primary 
  arma::vec observer_pos = this -> center_of_mass_vec[0] + this -> polydata_vec[0] -> GetLength() * 1E6 * observer_dir;
  arma::vec sun_pos = this -> center_of_mass_vec[0] + this -> polydata_vec[0] -> GetLength() * 1E6 * sun_dir;

  // Ray-tracing tolerance
  double tol = this -> polydata_vec[0] -> GetLength()/1E6;





  for (int body_index = 0; body_index < this -> number_of_bodies; ++body_index){

    vtkPolyData * input = this -> polydata_vec[body_index];

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds -> Allocate(VTK_CELL_SIZE);

    arma::vec::fixed<3> target_to_sun_dir_body_frame = BN_dcms_vec[body_index] * sun_dir;
    arma::vec::fixed<3> target_to_observer_dir_body_frame = BN_dcms_vec[body_index] * observer_dir;


  // The kept facets are then sampled and reverse ray-traced
    for (unsigned int facet_index = 0; facet_index != facets_in_view[body_index].size(); ++facet_index){

      double p0[3];
      double p1[3];
      double p2[3];

      input -> GetCellPoints(facets_in_view[body_index][facet_index],ptIds);
      input -> GetPoint(ptIds->GetId(0), p0);
      input -> GetPoint(ptIds->GetId(1), p1);
      input -> GetPoint(ptIds->GetId(2), p2);

      arma::vec::fixed<3> P0 = {p0[0],p0[1],p0[2]};
      arma::vec::fixed<3> P1 = {p1[0],p1[1],p1[2]};
      arma::vec::fixed<3> P2 = {p2[0],p2[1],p2[2]};

      arma::vec::fixed<3> n = arma::cross(P1 - P0, P2 - P0);

    // The number of points sampled from this facet is determined based on 
    // the relative size of this facet compared to the largest one in all the considered shapes
      int N_samples = int( N * arma::norm(n /2) / this -> min_area);

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


      // A random origin point is uniformly drawn from this facet
        arma::vec random = arma::randu<arma::vec>(2);
        double u = random(0);
        double v = random(1);
        arma::vec::fixed<3> origin = (1 - std::sqrt(u)) * P0 + std::sqrt(u) * ( 1 - v ) * P1 + std::sqrt(u) * v * P2 ;

    // Derived points are expressed in the body reference frame
      arma::vec point_above_surface = origin + 3 * tol * n; // the origin of the ray is moved 3*tol above the surface
      
      bool has_intersected = (this -> check_line_for_intersect(body_index,point_above_surface,sun_pos,BN_dcms_vec,positions_vec,tol)
        + this -> check_line_for_intersect(body_index,point_above_surface,observer_pos,BN_dcms_vec,positions_vec,tol));

      // If this point was not obscured, the return is weighed by the incidence on the inbound and outbout rays
      if (!has_intersected){
        measurements_temp[1] += cosi_sun * cosi_obs;
      }


    }
  }

}




}








