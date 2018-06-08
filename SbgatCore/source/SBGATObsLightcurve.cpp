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
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATObsLightcurve::~SBGATObsLightcurve(){
}

//----------------------------------------------------------------------------
// Description:
// This method computes the Cnm and Snm arrays of coefficients
// used in the spherical harmonics expansion of exterior gravity 
// about a constant-density polyhedron

int SBGATObsLightcurve::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){

  vtkInformation *inInfo =
  inputVector[0]->GetInformationObject(0);

  // call ExecuteData
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numCells, numPts;

  numCells = input->GetNumberOfCells();
  numPts = input->GetNumberOfPoints();
  if (numCells < 1 || numPts < 1){
    vtkErrorMacro( << "No data to measure...!");
    return 1;
  }

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds -> Allocate(VTK_CELL_SIZE);

  this -> bspTree =vtkSmartPointer<vtkModifiedBSPTree>::New();
  this -> bspTree->SetDataSet(input);

  // bspTree->SetMaxLevel(12);
  // bspTree->SetNumberOfCellsPerNode(16);
  this -> bspTree -> BuildLocator();

  // The center of mass of the input shape is saved
  vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();
  mass_filter -> SetInputData(input);
  mass_filter -> Update();
  this -> center_of_mass = mass_filter -> GetCenterOfMass();

  return 1;
}


void SBGATObsLightcurve::CollectMeasurementsSimpleSpin(
  std::vector<std::array<double, 2> > & measurements,
  const int & N,
  const double & dt,
  const double & period,
  const arma::vec & sun_dir,
  const arma::vec & observer_dir,
  const arma::vec & spin,
  const bool & penalize_indicence){

  // Containers
  std::vector<int> facets_in_view;
  vtkPolyData * input = vtkPolyData::SafeDownCast(this->GetInput(0));
  vtkIdType cellId, numCells, numPts, numIds;
  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds -> Allocate(VTK_CELL_SIZE);
  numCells = input -> GetNumberOfCells();
  numPts = input -> GetNumberOfPoints();


  // Angular velocity 
  double w = 2 * arma::datum::pi / period;
  arma::vec omega = spin * w;

  // Body-frame to inertial matrix
  arma::mat BN = RBK::prv_to_dcm(spin * dt * w );

  // Sun-to-target direction in body frame
  arma::vec target_to_sun_dir_body_frame = BN * arma::normalise(sun_dir);

  // Observer-to-target direction in body frame
  arma::vec target_to_observer_dir_body_frame = BN * arma::normalise(observer_dir);

  // Ray-tracing tolerance
  double tol = input -> GetLength()/1E6;


  // The surface area of the largest facet 
  double max_area = -1;
  // First, only facets that are in view of the sun
  // and the observer (based on their normal orientation) are kept
  for (cellId=0; cellId < numCells; cellId++){

    if ( input->GetCellType(cellId) != VTK_TRIANGLE){
      vtkWarningMacro(<< "Input data type must be VTK_TRIANGLE not "<< input->GetCellType(cellId));
      continue;
    }

    input -> GetCellPoints(cellId,ptIds);
    numIds = ptIds -> GetNumberOfIds();
    assert(numIds == 3);
    
    double p0[3];
    double p1[3];
    double p2[3];

    input->GetPoint(ptIds->GetId(0), p0);
    input->GetPoint(ptIds->GetId(1), p1);
    input->GetPoint(ptIds->GetId(2), p2);

    double e0[3];
    double e1[3];

    vtkMath::Subtract(p1,p0,e0);
    vtkMath::Subtract(p2,p0,e1);

    double n[3];
    
    vtkMath::Cross(e0,e1,n);
    max_area = std::max(max_area,vtkMath::Norm(n)/2);
    vtkMath::Normalize(n);

    // Check if the facet is in view of both the sun and the observer
    bool in_view = (vtkMath::Dot(n,target_to_sun_dir_body_frame.colptr(0)) > 0 && vtkMath::Dot(n,target_to_observer_dir_body_frame.colptr(0)) > 0);
    if (in_view){
      facets_in_view.push_back(cellId);
    }

  }



  // The vector holding the kept-facets is initialized
  std::array<double, 2>  measurements_temp;
  measurements_temp[0] = dt;
  measurements_temp[1] = 0;


  // The kept facets are then sampled and reverse ray-traced
  for (auto facet_index = 0; facet_index != facets_in_view.size(); ++facet_index){

    double p0[3];
    double p1[3];
    double p2[3];

    input -> GetCellPoints(facets_in_view[facet_index],ptIds);
    input -> GetPoint(ptIds->GetId(0), p0);
    input -> GetPoint(ptIds->GetId(1), p1);
    input -> GetPoint(ptIds->GetId(2), p2);

    arma::vec P0 = {p0[0],p0[1],p0[2]};
    arma::vec P1 = {p1[0],p1[1],p1[2]};
    arma::vec P2 = {p2[0],p2[1],p2[2]};

    arma::vec n = arma::cross(P1 - P0, P2 - P0);

     // The number of points sampled from this facet is determined based on 
  // the relative size of this facet compared to the largest one in the considered shape
    int N_facet = int( N * arma::norm(n /2) / max_area);

    // only need unit normal vector from here
    n = arma::normalise(arma::cross(P1 - P0, P2 - P0));

    double cosi_sun ;
    double cosi_obs ;

    if (penalize_indicence){
     cosi_sun = vtkMath::Dot(n.colptr(0),target_to_sun_dir_body_frame.colptr(0));
     cosi_obs = vtkMath::Dot(n.colptr(0),target_to_observer_dir_body_frame.colptr(0));

   }
   else{

    cosi_sun = 1;
    cosi_obs = 1;
  }


    // The incidence is constant across the facet

    // A maximum of N points are sampled from this facet
  // #pragma omp parallel for
  for (int i = 0; i < N_facet; ++i){

    bool keep_point = true;

      // A random impact point is uniformly drawn from this facet
    arma::vec random = arma::randu<arma::vec>(2);
    double u = random(0);
    double v = random(1);
    arma::vec impact =  (1 - std::sqrt(u)) * P0 + std::sqrt(u) * ( 1 - v ) * P1 + std::sqrt(u) * v * P2 ;





    // Checking if point is in view of sun
    {
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      vtkSmartPointer<vtkPoints> verts = vtkSmartPointer<vtkPoints>::New();

      arma::vec point_towards_sun = impact + input -> GetLength() * 1e3 * target_to_sun_dir_body_frame;

      this -> bspTree -> IntersectWithLine(impact.colptr(0), point_towards_sun.colptr(0), tol, verts, cellIds);

      // All the detected intersections are checked.
      for (int intersect_index = 0; intersect_index < verts -> GetNumberOfPoints(); ++intersect_index){

        double intersect[3];
        double intersect_to_origin[3];
        verts -> GetPoint(intersect_index,intersect);
        vtkMath::Subtract(intersect,point_towards_sun.colptr(0),intersect_to_origin);

        // If the intersect is between the radar and the impact point
        if (vtkMath::Norm(intersect_to_origin) + tol < arma::norm(impact - point_towards_sun)){

        // Reject this point
          keep_point = false;
          break;
        }

      }
    }


    // Checking if point is in view of the observer
    if(keep_point){
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      vtkSmartPointer<vtkPoints> verts = vtkSmartPointer<vtkPoints>::New();


      arma::vec point_towards_observer = impact + input -> GetLength() * 1e3 * target_to_observer_dir_body_frame;


      this -> bspTree -> IntersectWithLine(impact.colptr(0), point_towards_observer.colptr(0), tol, verts, cellIds);

      // All the detected intersections are checked.
      for (int intersect_index = 0; intersect_index < verts -> GetNumberOfPoints(); ++intersect_index){

        double intersect[3];
        double intersect_to_origin[3];
        verts -> GetPoint(intersect_index,intersect);
        vtkMath::Subtract(intersect,point_towards_observer.colptr(0),intersect_to_origin);

        // If the intersect is between the radar and the impact point
        if (vtkMath::Norm(intersect_to_origin) + tol < arma::norm(impact - point_towards_observer)){

        // Reject this point
          keep_point = false;
          break;
        }

      }
    }


    if (keep_point){

        // the return is weighed by the incidence on the inbound and outbout rays
      assert(cosi_sun * cosi_obs <= 1);
      assert(cosi_sun * cosi_obs >= 0);

      measurements_temp[1] += cosi_sun * cosi_obs;

    }

  }
}

measurements.push_back(measurements_temp);



}




void SBGATObsLightcurve::SaveLightCurveData(const std::vector<std::array<double, 2> > & measurements,
  std::string savepath){

  arma::mat time_luminosity_series(measurements.size(),2);

  for (int i = 0; i < measurements.size(); ++i){
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
