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

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SBGATShapeUncertainty.hpp"
#include "SBGATMassProperties.hpp"
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkCleanPolyData.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>
#include <json.hpp>
#include <ShapeUQLib/ShapeUQLib.hpp>


vtkStandardNewMacro(SBGATShapeUncertainty);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATShapeUncertainty::SBGATShapeUncertainty(){

  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATShapeUncertainty::~SBGATShapeUncertainty(){
}

//----------------------------------------------------------------------------
// Description:
// This method builds an internal duplicate of the provided shape model
// suitable for use in Sbgat's dependendy SBGATShapeUQLib

int SBGATShapeUncertainty::RequestData(vtkInformation* vtkNotUsed( request ),vtkInformationVector** inputVector,vtkInformationVector* vtkNotUsed( outputVector )){
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

 vtkPolyData * input_unclean = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkCleanPolyData> cleaner =
  vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner -> SetInputData (input_unclean);
  cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
  cleaner -> Update();

  vtkPolyData * input = cleaner -> GetOutput();

  vtkSmartPointer<SBGATMassProperties> center_of_mass_filter =
  vtkSmartPointer<SBGATMassProperties>::New();

  center_of_mass_filter -> SetInputData(input);
  center_of_mass_filter -> Update();

  arma::vec x = - center_of_mass_filter -> GetCenterOfMass();

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();

  transform -> Translate( x(0),  x(1),  x(2));

  vtkSmartPointer<vtkTransformPolyDataFilter> filter =
  vtkSmartPointer<vtkTransformPolyDataFilter>::New();

  filter -> SetInputData(input);
  filter -> SetTransform(transform);
  filter -> Update();

  input = filter -> GetOutput();


  vtkIdType numCells, numPts, numIds;

  numCells = input->GetNumberOfCells();
  numPts = input->GetNumberOfPoints();
  if (numCells < 1 || numPts < 1){
    vtkErrorMacro( << "No data to measure...!");
    return 1;
  }

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds -> Allocate(VTK_CELL_SIZE);

  std::vector<ControlPoint> vertices_coordinates(input -> GetNumberOfPoints());
  std::vector<std::vector<int > > facet_vertices_indices_coordinates(input -> GetNumberOfCells());
  std::vector<int> super_elements(input -> GetNumberOfCells());

  for(int i = 0; i < input -> GetNumberOfPoints(); ++i) {
    double coordinates[3];
    input -> GetPoint(i,coordinates);

    arma::vec::fixed<3> coords = {coordinates[0],coordinates[1],coordinates[2]};
    vertices_coordinates[i].set_point_coordinates(coordinates);

  } 


  for(int i = 0; i < numCells; ++i) {

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds -> Allocate(VTK_CELL_SIZE);

    input -> GetCellPoints(i,ptIds);

    facet_vertices_indices_coordinates[i].push_back(ptIds -> GetId(0));
    facet_vertices_indices_coordinates[i].push_back(ptIds -> GetId(1));
    facet_vertices_indices_coordinates[i].push_back(ptIds -> GetId(2));

    super_elements[i] = i;
  }


  ShapeModelTri<ControlPoint> shape_model_tri(facet_vertices_indices_coordinates,
    super_elements,vertices_coordinates);

  this -> shape_model_bezier = std::make_shared<ShapeModelBezier<ControlPoint>>(ShapeModelBezier<ControlPoint>(ShapeModelBezier<ControlPoint>(shape_model_tri,"",nullptr)));




  for (int e = 0; e < this -> shape_model_bezier -> get_NElements(); ++e){

    this -> shape_model_bezier -> get_element(e).set_owning_shape(this -> shape_model_bezier.get());

  }

  return 1;
}


void SBGATShapeUncertainty::ComputeInertiaStatistics(double l,double sigma){
  this -> Update();

  this -> shape_model_bezier -> compute_point_covariances(std::pow(sigma,2),l);

  this -> shape_model_bezier -> compute_all_statistics();

  this -> volume_variance = std::pow(this -> shape_model_bezier -> get_volume_sd(),2);
  this -> center_of_mass_covariance = this -> shape_model_bezier -> get_cm_cov();
  this -> inertia_tensor_parametrization_covariance = this -> shape_model_bezier -> get_inertia_cov();
  this -> principal_dimensions_covariance = this -> shape_model_bezier -> get_P_dims();
  this -> principal_moments_covariance = this -> shape_model_bezier -> get_P_moments();
  this -> mrp_covariance = this -> shape_model_bezier ->  get_mrp_cov();

}

void SBGATShapeUncertainty::ComputeInertiaStatisticsMC(int N_samples,double l,double sigma,std::string output_dir){
  
  this -> Update();
  this -> shape_model_bezier -> compute_point_covariances(std::pow(sigma,2),l);
  this -> shape_model_bezier -> compute_shape_covariance_sqrt();

  arma::vec results_volume;
  arma::mat results_cm,results_inertia,results_moments,results_mrp,results_dims;

  this -> shape_model_bezier -> run_monte_carlo(N_samples,
    results_volume,
    results_cm,
    results_inertia,
    results_moments,
    results_mrp,
    results_dims,
    output_dir);


  arma::vec results_cm_mean = arma::mean(results_cm,1);
  arma::vec results_inertia_mean = arma::mean(results_inertia,1);
  arma::vec results_moments_mean = arma::mean(results_moments,1);
  arma::vec results_dims_mean = arma::mean(results_dims,1);
  arma::vec results_mrp_mean = arma::mean(results_mrp,1);

  arma::mat cov_cm_mc = arma::zeros(3,3);
  arma::mat cov_inertia_mc = arma::zeros(6,6);
  arma::mat cov_moments_mc = arma::zeros(4,4);
  arma::mat cov_dims_mc = arma::zeros<arma::mat>(3,3);
  arma::mat cov_mrp_mc = arma::zeros(3,3);
  
  for (unsigned int i = 0; i < results_cm.n_cols; ++i){

    cov_cm_mc +=  (results_cm.col(i) - results_cm_mean) * (results_cm.col(i) - results_cm_mean).t();
    cov_inertia_mc +=  (results_inertia.col(i) - results_inertia_mean) * (results_inertia.col(i) - results_inertia_mean).t();
    cov_moments_mc +=  (results_moments.col(i) - results_moments_mean) * (results_moments.col(i) - results_moments_mean).t();
    cov_dims_mc +=  (results_dims.col(i) - results_dims_mean) * (results_dims.col(i) - results_dims_mean).t();
    cov_mrp_mc +=  (RBK::dcm_to_mrp(RBK::mrp_to_dcm(results_mrp.col(i))*RBK::mrp_to_dcm(results_mrp_mean).t())) * (RBK::dcm_to_mrp(RBK::mrp_to_dcm(results_mrp.col(i))*RBK::mrp_to_dcm(results_mrp_mean).t())).t();

  }

  cov_cm_mc *= 1./(results_cm.n_cols-1);
  cov_inertia_mc *= 1./(results_inertia.n_cols-1);
  cov_moments_mc *= 1./(results_moments.n_cols-1);
  cov_dims_mc *= 1./(results_dims.n_cols-1);
  cov_mrp_mc *= 1./(results_mrp.n_cols-1);

  this -> volume_variance = std::pow(arma::stddev(results_volume,1),2);
  this -> center_of_mass_covariance = cov_cm_mc;
  this -> inertia_tensor_parametrization_covariance = cov_inertia_mc;
  this -> principal_dimensions_covariance = cov_dims_mc;
  this -> principal_moments_covariance = cov_moments_mc;

  this -> mrp_covariance = cov_mrp_mc;



}



void SBGATShapeUncertainty::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATShapeUncertainty::PrintTrailer(ostream& os, vtkIndent indent) {

}


//----------------------------------------------------------------------------
void SBGATShapeUncertainty::PrintSelf(std::ostream& os, vtkIndent indent){

  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
  if (!input){
    return;
  }

}





