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
  Module:    SBGATObsRadar.cpp

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include <SBGATObsRadar.hpp>
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



vtkStandardNewMacro(SBGATObsRadar);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATObsRadar::SBGATObsRadar(){

  this -> SetNumberOfOutputPorts(0);
  this -> SetNumberOfInputPorts(2);

  this -> SetInputData(0,nullptr);
  this -> SetInputData(1,nullptr);



}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATObsRadar::~SBGATObsRadar(){
}






int SBGATObsRadar::FillInputPortInformation( int port, vtkInformation* info ){
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




int SBGATObsRadar::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){

  this -> bspTree_vec.clear();
  this -> polydata_vec.clear();

  vtkSmartPointer<SBGATMassProperties> mass_filter = vtkSmartPointer<SBGATMassProperties>::New();


  // Processing the primary
  vtkInformation *inInfo0 = inputVector[0]->GetInformationObject(0);
  vtkPolyData * primary = vtkPolyData::SafeDownCast(inInfo0->Get(vtkDataObject::DATA_OBJECT()));

  

  vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
  tree -> SetDataSet(primary);
  tree -> BuildLocator();
  
  this -> bspTree_vec.push_back(tree);
  this -> polydata_vec.push_back(primary);

  mass_filter -> SetInputData(primary);
  mass_filter -> Update();
  this -> center_of_mass_vec.push_back(mass_filter -> GetCenterOfMass());


  // Processing the secondary, if any
  vtkInformation *inInfo1 = inputVector[1]->GetInformationObject(0);

  if(inInfo1->Get(vtkDataObject::DATA_OBJECT()) != nullptr){
   vtkPolyData * secondary =  vtkPolyData::SafeDownCast(inInfo1->Get(vtkDataObject::DATA_OBJECT()));
   vtkSmartPointer<vtkModifiedBSPTree> tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
   tree -> SetDataSet(secondary);
   tree -> BuildLocator();
   this -> bspTree_vec.push_back(tree);
   this -> polydata_vec.push_back(secondary);

   mass_filter -> SetInputData(secondary);
   mass_filter -> Update();
   this -> center_of_mass_vec.push_back(mass_filter -> GetCenterOfMass());

 }

 this -> number_of_bodies = polydata_vec.size();
 
 // The surface area of the largest facet amongst all considered shapes is found
 this -> find_max_facet_surface_area();

 return 1;
}



void SBGATObsRadar::CollectMeasurementsSimpleSpin(SBGATRadarObsSequence & measurements_sequence,  
  const int & N,
  const double & dt,
  const std::vector<double> & period_vec,
  const arma::vec & radar_dir,
  const std::vector<arma::vec> & positions_vec,
  const std::vector<arma::vec> & velocities_vec,
  const std::vector<arma::vec> & spin_vec,
  const bool & penalize_indicence){


  if(period_vec.size() != spin_vec.size() || period_vec.size() != positions_vec.size() || period_vec.size() != velocities_vec.size())
    throw(std::runtime_error("Incompatible input dimensions"));


  // Containers
  std::vector<std::vector<int> > facets_in_view;
  std::vector<arma::mat> BN_dcms_vec ;
  std::vector<arma::vec> omega_vec ;


  // Pre-allocating for all inputs  
  for (int i = 0; i < this -> number_of_bodies; ++i){
    double w = 2 * arma::datum::pi / period_vec[i];
    facets_in_view.push_back(std::vector<int>());
    BN_dcms_vec.push_back(RBK::prv_to_dcm(spin_vec[i] * dt * w ));
    omega_vec.push_back(spin_vec[i] * w);
  }


  // First, only facets that are in view of the radar (based on their normal orientation) are kept
  // Then, the facets that are potentially in view are ray-traced to the observer,
  // checking with potential intersects with all bodies

  std::vector<arma::vec> dir_to_check_vec = {radar_dir};

  for (int i = 0; i < this -> number_of_bodies; ++i){
    this -> prefind_facets_inview(facets_in_view[i],i,dir_to_check_vec,BN_dcms_vec,positions_vec);
  }


  // The vector holding the kept-facets is initialized
  std::vector<std::array<double, 3> > measurements_temp;
  
  this -> reverse_ray_trace(measurements_sequence,facets_in_view,radar_dir,N,penalize_indicence,
    BN_dcms_vec,positions_vec,velocities_vec,omega_vec);



}


void SBGATObsRadar::reverse_ray_trace(SBGATRadarObsSequence & measurements_sequence,
  const std::vector<std::vector<int> > & facets_in_view,
  const arma::vec & radar_dir,
  const int N,
  const bool penalize_indicence,
  const std::vector<arma::mat> & BN_dcms_vec,
  const std::vector<arma::vec> & positions_vec,
  const std::vector<arma::vec> & velocities_vec,
  const std::vector<arma::vec> & omega_vec){

  std::vector<std::array<double, 3> > measurements;

  // Ray-tracing tolerance
  double tol = this -> polydata_vec[0] -> GetLength()/1E6;

  // The radar is positionned with respect to the primary 
  arma::vec radar_pos = this -> center_of_mass_vec[0] + this -> polydata_vec[0] -> GetLength() * 1E6 * radar_dir;

  for (int body_index = 0; body_index < this -> number_of_bodies; ++body_index){

    vtkPolyData * input = this -> polydata_vec[body_index];

    vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
    ptIds -> Allocate(VTK_CELL_SIZE);

    arma::vec target_to_radar_dir_body_frame = BN_dcms_vec[body_index] * radar_dir;




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
      int N_samples = int( N * arma::norm(n /2) / this -> min_area);

    // only need unit normal vector from here
      n = arma::normalise(arma::cross(P1 - P0, P2 - P0));

      double cosi_radar ;

    // Computing the ray incidence at impact if needed
      if (penalize_indicence){
        cosi_radar = vtkMath::Dot(n.colptr(0),target_to_radar_dir_body_frame.colptr(0));
      }
      else{
        cosi_radar = 1;
      }

    // A maximum of N points are sampled from this facet
    // #pragma omp parallel for
      for (int i = 0; i < N_samples; ++i){

      // A random origin point is uniformly drawn from this facet
        arma::vec random = arma::randu<arma::vec>(2);
        double u = random(0);
        double v = random(1);
        arma::vec origin = (1 - std::sqrt(u)) * P0 + std::sqrt(u) * ( 1 - v ) * P1 + std::sqrt(u) * v * P2 ;

      arma::vec point_above_surface = origin + 3 * tol * n; // the origin of the ray is moved 3*tol above the surface

      bool has_intersected = this -> check_line_for_intersect(body_index,point_above_surface,radar_pos,BN_dcms_vec,positions_vec,tol);

      // If this point was not obscured, the return is weighed by the incidence on the inbound and outbout rays

      if (!has_intersected){

        // origin is the sampled point, expressed in the body reference frame
        // Its actual position should be in the inertial frame, accounting for the assigned position of the body

        arma::vec origin_cm = BN_dcms_vec[body_index].t() * (origin - this -> center_of_mass_vec[body_index]);
        arma::vec origin_inertial = origin_cm + positions_vec[body_index];

        // Range
        double range = arma::norm(origin_inertial - radar_pos);

        // The velocity at the impact point is a combination of the orbital and rotational velocities
        arma::vec velocity = velocities_vec[body_index] + arma::cross(omega_vec[body_index],origin_cm);
        
        double range_rate = arma::dot(origin_inertial - radar_pos,velocity) / range;

        std::array<double, 3> measurement = {{range,range_rate,std::pow(cosi_radar,2)}};
        measurements.push_back(measurement);
      }


    }
  }
}

measurements_sequence.push_back(measurements);

}

void SBGATObsRadar::BinObservations(
  const SBGATRadarObsSequence & measurements_sequence,
  const double & r_bin,
  const double & rr_bin){

  // The container holding the images is pre-allocated
  this -> images.clear();
  this -> max_value = -1;

  for (unsigned int i = 0; i < measurements_sequence.size(); ++i){

    auto measurements = measurements_sequence[i];

    this -> images.push_back(vtkSmartPointer<vtkImageData>::New());

  // Using arma::vec to get statistics
    arma::vec ranges(measurements.size());
    arma::vec range_rates(measurements.size());

    for (unsigned int i = 0; i < measurements.size(); ++i){
      ranges(i) = measurements[i][0];
      range_rates(i) = measurements[i][1];
    }

  // Subtracting the mean range
    ranges -= arma::mean(ranges);

  // The extent of the data is extracted
    double r_extent = std::abs(ranges.max() - ranges.min()) * this -> scaleFactor;
    double rr_extent = std::abs(range_rates.max() - range_rates.min()) * this -> scaleFactor;


  // The bin counts are determined from the specified bin sizes and data extent
    int n_bin_r = (int)(r_extent / r_bin);
    int n_bin_rr = (int)(rr_extent / rr_bin);

    // Checking if one dimension is "empty" (i.e has zero bins)
    if (n_bin_r == 0 || n_bin_rr == 0){
      throw(std::runtime_error("The prescribed bin sizes yielded " + std::to_string(n_bin_r) + " range bins and " + std::to_string(n_bin_rr) + " range-rate bins. "));
    }

    // The image is formed by "binning in" the measurements

    this -> images[i] -> SetExtent(0, n_bin_rr - 1, 0, n_bin_r - 1, 0, 0);
    this -> images[i] -> AllocateScalars(VTK_DOUBLE, 1);

    double *dPtr = static_cast<double *>(this -> images[i]->GetScalarPointer(0, 0, 0));

    // The image is initialized 
    #pragma omp parallel for
    for (int row = 0; row < n_bin_r; ++row){
      for (int col = 0; col < n_bin_rr; ++col){
        int k = col + row * n_bin_rr;
        dPtr[k] = 0;
      }
    }


    double max_range = ranges.max();
    double max_range_rate = range_rates.max();

    // The histogram is built
    #pragma omp parallel for
    for (unsigned int mes = 0; mes < measurements.size(); ++mes){

      // Flipping the image
      int row = int( ( - ranges(mes) +  max_range)/ r_bin );
      int col = int( ( - range_rates(mes) + max_range_rate )/ rr_bin );

      // The last bin is inclusive
      if (row == n_bin_r){
        --row;
      }
      if (col == n_bin_rr){
        --col;
      }

      int k = col + row * n_bin_rr;

      // Adding the cosine of the incidence 
      // If the incidence is 0 then measurements[mes][2] == 1. 
      // if the incidence is 90 deg then measurements[mes][2] == 0
      // if the incidence is not used then measurements[mes][2] == 1 always
      dPtr[k] += measurements[mes][2];

    }

    // The maximum image value is extracted
    vtkDataArray * scalars = this -> images[i] -> GetPointData() -> GetScalars();

    for (vtkIdType tupleIdx = 0; tupleIdx < scalars -> GetNumberOfTuples(); ++tupleIdx){
     this -> max_value = std::max(max_value,scalars -> GetTuple1(tupleIdx));
   }

 }

}

void SBGATObsRadar::SaveImages( std::string savepath){

  vtkSmartPointer<vtkImageCast> cast = vtkSmartPointer<vtkImageCast>::New();
  vtkSmartPointer<vtkPNGWriter> writer = vtkSmartPointer<vtkPNGWriter> ::New();

  for (unsigned int i = 0; i < this -> images.size(); ++i){

   // Each image is normalized before being saved
   vtkSmartPointer<vtkImageData> image = vtkSmartPointer<vtkImageData>::New();
   image -> DeepCopy(this -> images[i]);

   vtkDataArray * scalars = image -> GetPointData() -> GetScalars();

   for (vtkIdType tupleIdx = 0; tupleIdx < scalars -> GetNumberOfTuples(); ++tupleIdx){
    scalars -> SetTuple1(tupleIdx,scalars -> GetTuple1(tupleIdx) * 255 / this -> max_value );
  }


  cast -> SetInputData(image);
  cast -> SetOutputScalarTypeToUnsignedChar();
  std::string complete_save_path = savepath + "image_" + std::to_string(i) + ".png";

  writer -> SetFileName(complete_save_path.c_str());
  writer -> SetInputConnection(cast->GetOutputPort());
  writer -> Write();
}

}

std::vector<vtkSmartPointer<vtkImageData>> SBGATObsRadar::GetImages() const{
  return this -> images;

}

void SBGATObsRadar::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATObsRadar::PrintTrailer(ostream& os, vtkIndent indent) {

}


//----------------------------------------------------------------------------
void SBGATObsRadar::PrintSelf(std::ostream& os, vtkIndent indent){

  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
  if (!input)
  {
    return;
  }

}
