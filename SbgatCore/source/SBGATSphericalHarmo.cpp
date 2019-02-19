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
#include <SBGATSphericalHarmo.hpp>
#include <SBGATMassProperties.hpp>
#include <SHARMLib.hpp>
#include <json.hpp>
#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkCleanPolyData.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>
vtkStandardNewMacro(SBGATSphericalHarmo);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATSphericalHarmo::SBGATSphericalHarmo(){

  this -> normalized = false;
  this -> degreeSet = false;
  this -> densitySet = false;
  this -> referenceRadiusSet = false;
  this -> scaleFactorSet = false;
  this -> setFromJSON = false;


  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATSphericalHarmo::~SBGATSphericalHarmo(){
}

//----------------------------------------------------------------------------
// Description:
// This method computes the Cnm and Snm arrays of coefficients
// used in the spherical harmonics expansion of exterior gravity 
// about a constant-density polyhedron

int SBGATSphericalHarmo::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){

  if (this -> setFromJSON){
    return 1;
  }


  vtkInformation *inInfo =
  inputVector[0]->GetInformationObject(0);

  // call ExecuteData
  vtkPolyData * input_unclean = vtkPolyData::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkSmartPointer<vtkCleanPolyData> cleaner =
  vtkSmartPointer<vtkCleanPolyData>::New();
  cleaner -> SetInputData (input_unclean);
  cleaner -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );
  cleaner -> Update();

  vtkPolyData * input = cleaner -> GetOutput();
  vtkIdType cellId, numCells, numPts, numIds;

  numCells = input->GetNumberOfCells();
  numPts = input->GetNumberOfPoints();
  if (numCells < 1 || numPts < 1)
  {
    vtkErrorMacro( << "No data to measure...!");
    return 1;
  }

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds->Allocate(VTK_CELL_SIZE);


  if (!(this -> degreeSet && this -> densitySet && this -> referenceRadiusSet && this -> scaleFactorSet)){
    throw(std::runtime_error("Trying to evaluate spherical harmonics although the degree, density, reference radius and scale factor may have not been properly set"));
  }
  
  // Initialize variables 
  this -> Cnm.clear();
  this -> Snm.clear();
  this -> n_facets = numCells;
  this -> n_vertices = numPts;

  this -> Cnm = arma::zeros<arma::mat>(this -> degree + 1  , this -> degree + 1);
  this -> Snm = arma::zeros<arma::mat>(this -> degree + 1  , this -> degree + 1);
  
  vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
  mass_properties -> SetInputData(input);
  mass_properties -> Update();

  // Check that the shape is topologically closed
  assert(mass_properties -> CheckClosed());

  this -> totalMass = mass_properties -> GetVolume() * this -> density * std::pow(this -> scaleFactor,3);

  // Looping over all facets
  for (cellId=0; cellId < numCells; cellId++){

    if ( input->GetCellType(cellId) != VTK_TRIANGLE){
      vtkWarningMacro(<< "Input data type must be VTK_TRIANGLE not " << input->GetCellType(cellId));
      continue;
    }

    input->GetCellPoints(cellId,ptIds);
    numIds = ptIds->GetNumberOfIds();
    assert(numIds == 3);

    double r0[3];
    double r1[3];
    double r2[3];

    input->GetPoint(ptIds->GetId(0), r0);
    input->GetPoint(ptIds->GetId(1), r1);
    input->GetPoint(ptIds->GetId(2), r2);

    arma::mat Cnm2f(this -> degree + 1, this -> degree + 1);
    arma::mat Snm2f(this -> degree + 1, this -> degree + 1);

    // Call to SHARMLib here
    SHARMLib::ComputePolyhedralCS(
      Cnm2f,
      Snm2f,
      this -> degree,
      this -> referenceRadius,
      &r0[0],
      &r1[0],
      &r2[0],
      this -> normalized
      );

    this -> Cnm += Cnm2f ;
    this -> Snm += Snm2f ;

  }

  this -> Cnm /= mass_properties -> GetVolume();
  this -> Snm /= mass_properties -> GetVolume();

  return 1;
}








arma::vec::fixed<3> SBGATSphericalHarmo::GetAcceleration(const arma::vec::fixed<3> & pos){

  try{

    this -> Update();

    int n_max = 50;

    arma::mat b_bar_real = arma::zeros<arma::mat>(n_max + 3,n_max + 3);
    arma::mat b_bar_imag = arma::zeros<arma::mat>(n_max + 3,n_max + 3);

    SHARMLib::GetBnmNormalizedExterior(this -> degree,
      b_bar_real,
      b_bar_imag,
      pos,
      this -> referenceRadius);

    double mu = this -> totalMass * arma::datum::G;

    double K0 = 0.5 * mu / std::pow(this -> referenceRadius * this -> scaleFactor,2);
    double x_ddot = 0;
    double y_ddot = 0;
    double z_ddot = 0;

    for (unsigned int nn = 0; nn <= this -> degree; nn++){

      double n = (double) nn;

      for (unsigned int mm = 0; mm<=nn; mm++){

        double m = (double) mm;
        double delta_1_m;

        if (mm == 1){
          delta_1_m = 1.0;
        }
        else{
          delta_1_m = 0.0;
        } 

        double K1 = sqrt( (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+3.0) );
        double K2 = sqrt( (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) );
        double K3 = sqrt( 2.0 * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+3.0) );

        if (mm == 0){

          x_ddot -= 2.0 * K0 * ( this -> Cnm(nn,mm) * K1 * b_bar_real(nn+1,mm+1) );
          y_ddot -= 2.0 * K0 * ( this -> Cnm(nn,mm) * K1 * b_bar_imag(nn+1,mm+1) );
          z_ddot -= 2.0 * K0 * ( this -> Cnm(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0)) * b_bar_real(nn+1,mm) );

        }
        else{

          x_ddot += K0 * ( -this -> Cnm(nn,mm) * K2 * b_bar_real(nn+1,mm+1) -this -> Snm(nn,mm) * K2 * b_bar_imag(nn+1,mm+1) +this -> Cnm(nn,mm) * K3 * b_bar_real(nn+1,mm-1) +this -> Snm(nn,mm) * K3 * b_bar_imag(nn+1,mm-1));
          y_ddot += K0 * ( -this -> Cnm(nn,mm) * K2 * b_bar_imag(nn+1,mm+1) +this -> Snm(nn,mm) * K2 * b_bar_real(nn+1,mm+1) -this -> Cnm(nn,mm) * K3 * b_bar_imag(nn+1,mm-1) +this -> Snm(nn,mm) * K3 * b_bar_real(nn+1,mm-1));
          z_ddot -= 2.0 * K0 * ( this -> Cnm(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2.0*n+1.0)/(2.0*n+3.0)) * b_bar_real(nn+1,mm) +this -> Snm(nn,mm)*sqrt((n-m+1.0)*(n+m+1.0)*(2*n+1.0)/(2.0*n+3.0)) * b_bar_imag(nn+1,mm) );

        } 

      } 

    } 

    arma::vec::fixed<3> acceleration = {x_ddot,y_ddot,z_ddot};

    acceleration *= 1./this -> scaleFactor;

    return acceleration;

  }

  catch(std::runtime_error & error){

    std::cout << "an std::runtime_error occured inside SBGATSphericalHarmo::GetAcceleration. returning (0,0,0)\n";
    return arma::zeros<arma::vec>(3);

  }

} 



void SBGATSphericalHarmo::GetGravityGradientMatrix(const arma::vec::fixed<3> & pos,
  arma::mat::fixed<3,3> & dAccdPos){

  try{

    this -> Update();

    int n_max = 50;

    arma::mat b_bar_real = arma::zeros<arma::mat>(n_max + 3,n_max + 3);
    arma::mat b_bar_imag = arma::zeros<arma::mat>(n_max + 3,n_max + 3);

    SHARMLib::GetBnmNormalizedExterior(this -> degree,
      b_bar_real,
      b_bar_imag,
      pos,
      this -> referenceRadius);

    double mu = this -> totalMass * arma::datum::G;

    double K0 = 0.25 * mu / std::pow(this -> referenceRadius * this -> scaleFactor,3);

    double ddU_dxdx = 0;
    double ddU_dydy = 0;
    double ddU_dxdy = 0;
    double ddU_dxdz = 0;
    double ddU_dydz = 0;
    double ddU_dzdz = 0;

    for (unsigned int nn = 0; nn <= this -> degree; nn++){

      double n = (double) nn;

      for (unsigned int mm = 0; mm<=nn; mm++){

        double m = (double) mm;
        double delta_1_m,delta_2_m;

        if (mm == 1){
          delta_1_m = 1.0;
        }
        else{
          delta_1_m = 0.0;
        } 

        if (mm == 2){
          delta_2_m = 1.0;
        }
        else{
          delta_2_m = 0.0;
        }

        double K1 = sqrt( (n+m+4.0) * (n+m+3.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
        double K2 = sqrt( (n-m+2.0) * (n-m+1.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
        double K3 = sqrt( 2.0 * (n-m+4.0) * (n-m+3.0) * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_2_m) / (2.0*n+5.0) );
        double K4 = sqrt( (n+5.0) * (n+4.0) * (n+3.0) * (n+2.0) * (2.0*n+1.0) / (2.0*n+5.0) );
        double K5 = sqrt( (n+3.0) * (n+2.0) * (n+1.0) * n * (2.0*n+1.0) / (2.0*n+5.0) );
        double K6 = sqrt( (n+4.0) * (n+3.0) * (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+5.0) );
        double K7 = sqrt( (2.0*n+1.0) / (2.0*n+5.0) );
        double K8 = sqrt( (n-m+1.0) * (n+m+3.0) * (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+5.0) );
        double K9 = sqrt( 2.0 * (n+m+1.0) * (n-m+3.0) * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+5.0) );
        double K10= sqrt( (n+3.0) * (n+2.0) * (n+1.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+5.0) );


        /*// Partial expressions */
        if (mm == 0){

          ddU_dxdx += 2.0 * K0 * ( this -> Cnm(nn,mm) * K6 * b_bar_real(nn+2,mm+2) -(n+2)*(n+1) * K7*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) );
          ddU_dydy -= 2.0 * K0 * ( this -> Cnm(nn,mm) * K6 * b_bar_real(nn+2,mm+2) +(n+2)*(n+1) * K7*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) );
          ddU_dxdy += 2.0 * K0 * ( this -> Cnm(nn,mm) * K6 * b_bar_imag(nn+2,mm+2) );
          ddU_dxdz += 4.0 * K0 * ( this -> Cnm(nn,mm) * K10 * b_bar_real(nn+2,mm+1) );
          ddU_dydz += 4.0 * K0 * ( this -> Cnm(nn,mm) * K10 * b_bar_imag(nn+2,mm+1) );
          ddU_dzdz += 4.0 * K0 * ( this -> Cnm(nn,mm) * K2  * b_bar_real(nn+2,mm) );

        }
        else if (mm == 1){

          ddU_dxdx += K0 * ( this -> Cnm(nn,mm) * K4 * b_bar_real(nn+2,mm+2) +this -> Snm(nn,mm) * K4 * b_bar_imag(nn+2,mm+2) -3.0 * K5*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) -  K5*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm));
          ddU_dydy -= K0 * ( this -> Cnm(nn,mm) * K4 * b_bar_real(nn+2,mm+2) +this -> Snm(nn,mm) * K4 * b_bar_imag(nn+2,mm+2) +  K5*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) +3.0 * K5*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm));
          ddU_dxdy -= K0 * ( this -> Snm(nn,mm) * K4 * b_bar_real(nn+2,mm+2) -this -> Cnm(nn,mm) * K4 * b_bar_imag(nn+2,mm+2) +  K5*this -> Snm(nn,mm) * b_bar_real(nn+2,mm) +  K5*this -> Cnm(nn,mm) * b_bar_imag(nn+2,mm));
          ddU_dxdz += 2.0 * K0 * ( this -> Cnm(nn,mm) * K8 * b_bar_real(nn+2,mm+1) +this -> Snm(nn,mm) * K8 * b_bar_imag(nn+2,mm+1) -  K9*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm-1) -K9*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm-1) );
          ddU_dydz -= 2.0 * K0 * ( this -> Snm(nn,mm) * K8 * b_bar_real(nn+2,mm+1) -this -> Cnm(nn,mm) * K8 * b_bar_imag(nn+2,mm+1) +  K9*this -> Snm(nn,mm) * b_bar_real(nn+2,mm-1) -K9*this -> Cnm(nn,mm) * b_bar_imag(nn+2,mm-1) );
          ddU_dzdz += 4.0 * K0 * ( this -> Cnm(nn,mm) * K2 * b_bar_real(nn+2,mm)   +this -> Snm(nn,mm) * K2 * b_bar_imag(nn+2,mm));
        }
        else{

          ddU_dxdx += K0 * (  this -> Cnm(nn,mm) * K1 * b_bar_real(nn+2,mm+2) +this -> Snm(nn,mm) * K1 * b_bar_imag(nn+2,mm+2) -2.0 * K2*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) -2.0 * K2*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm) +K3*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm-2) +K3*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm-2));
          ddU_dydy -= K0 * (  this -> Cnm(nn,mm) * K1 * b_bar_real(nn+2,mm+2) +this -> Snm(nn,mm) * K1 * b_bar_imag(nn+2,mm+2) +2.0 * K2*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm) +2.0 * K2*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm) +K3*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm-2) +K3*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm-2));
          ddU_dxdy += K0 * ( -this -> Snm(nn,mm) * K1 * b_bar_real(nn+2,mm+2) +this -> Cnm(nn,mm) * K1 * b_bar_imag(nn+2,mm+2) +K3*this -> Snm(nn,mm) * b_bar_real(nn+2,mm-2) -K3*this -> Cnm(nn,mm) * b_bar_imag(nn+2,mm-2));
          ddU_dxdz += 2.0 * K0 * ( this -> Cnm(nn,mm) * K8 * b_bar_real(nn+2,mm+1) +this -> Snm(nn,mm) * K8 * b_bar_real(nn+2,mm+1) -  K9*this -> Cnm(nn,mm) * b_bar_real(nn+2,mm-1) -K9*this -> Snm(nn,mm) * b_bar_imag(nn+2,mm-1) );
          ddU_dydz -= 2.0 * K0 * ( this -> Snm(nn,mm) * K8 * b_bar_real(nn+2,mm+1) -this -> Cnm(nn,mm) * K8 * b_bar_imag(nn+2,mm+1) +  K9*this -> Snm(nn,mm) * b_bar_real(nn+2,mm-1) -K9*this -> Cnm(nn,mm) * b_bar_imag(nn+2,mm-1) );
          ddU_dzdz += 4.0 * K0 * ( this -> Cnm(nn,mm) * K2 * b_bar_real(nn+2,mm)   +this -> Snm(nn,mm) * K2 * b_bar_imag(nn+2,mm));

        }

      } 

    } 

    dAccdPos(0,0) = ddU_dxdx;
    dAccdPos(1,1) = ddU_dydy;
    dAccdPos(2,2) = ddU_dzdz;

    dAccdPos(0,1) = ddU_dxdy;
    dAccdPos(1,0) = ddU_dxdy;

    dAccdPos(0,2) = ddU_dxdz;
    dAccdPos(2,0) = ddU_dxdz;

    dAccdPos(1,2) = ddU_dydz;
    dAccdPos(2,1) = ddU_dydz;


  }

  catch(std::runtime_error & error){

    std::cout << "an std::runtime_error occured inside SBGATSphericalHarmo::GetGravityGradientMatrix. returning zero matrix\n";
    dAccdPos = arma::zeros<arma::mat>(3,3);

  }

} 

void SBGATSphericalHarmo::GetPartialHarmonics(const arma::vec::fixed<3> & pos,
  arma::mat & partial_C, 
  arma::mat & partial_S){

  int Ccounter = 0;
  int Scounter = 0;

  int n_max = 50;
  double mu = this -> totalMass * arma::datum::G;

  double K0 = 0.5 * mu / std::pow(this -> referenceRadius * this -> scaleFactor,2);
  
  arma::mat b_bar_real = arma::zeros<arma::mat>(n_max + 3,n_max + 3);
  arma::mat b_bar_imag = arma::zeros<arma::mat>(n_max + 3,n_max + 3);

  SHARMLib::GetBnmNormalizedExterior(this -> degree,
    b_bar_real,
    b_bar_imag,
    pos,
    this -> referenceRadius);


  partial_C.set_size(3,static_cast<int>((this -> degree + 1) * (this -> degree + 2)/2));
  partial_S.set_size(3,static_cast<int>((this -> degree + 1) * (this -> degree + 2)/2 - static_cast<int>(this -> degree + 1)));

  for (unsigned int nn = 0; nn <= this -> degree; nn++){

    double n = (double) nn;

    for (unsigned int mm = 0; mm<=nn; mm++){

      double m = (double) mm;
      double delta_1_m;

      if (mm == 1){
        delta_1_m = 1.0;
      }
      else{
        delta_1_m = 0.0;
      } 

      double K1 = sqrt( (n+2.0) * (n+1.0) * (2.0*n+1.0) / 2.0 / (2.0*n+3.0) );
      double K2 = sqrt( (n+m+2.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) );
      double K3 = sqrt( 2.0 * (n-m+2.0) * (n-m+1.0) * (2.0*n+1.0) / (2.0 - delta_1_m) / (2.0*n+3.0) );

      if (mm == 0){

        partial_C(0,Ccounter)   = - 2.0 * K0 * ( K1 * b_bar_real(nn+1,mm+1) );
        partial_C(1,Ccounter) = - 2.0 * K0 * ( K1 * b_bar_imag(nn+1,mm+1) );
        partial_C(2,Ccounter) = - 2.0 * K0 * ( sqrt( (n-m+1.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) ) * b_bar_real(nn+1,mm) );
        Ccounter += 1;

      }           
      else{

        partial_C(0,Ccounter)   = K0 * ( -K2 * b_bar_real(nn+1,mm+1) +K3 * b_bar_real(nn+1,mm-1) );
        partial_C(1,Ccounter) = K0 * ( -K2 * b_bar_imag(nn+1,mm+1) -K3 * b_bar_imag(nn+1,mm-1) );
        partial_C(2,Ccounter) = -2.0 * K0 * ( sqrt( (n-m+1.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) ) * b_bar_real(nn+1,mm) );
        Ccounter += 1;
        
        partial_S(0,Scounter)   = K0 * ( -K2 * b_bar_imag(nn+1,mm+1) + K3 * b_bar_imag(nn+1,mm-1) );
        partial_S(1,Scounter) = K0 * (  K2 * b_bar_real(nn+1,mm+1) + K3 * b_bar_real(nn+1,mm-1) );
        partial_S(2,Scounter) = -2.0 * K0 * ( sqrt( (n-m+1.0) * (n+m+1.0) * (2.0*n+1.0) / (2.0*n+3.0) ) * b_bar_imag(nn+1,mm) );
        Scounter += 1;

      } 

    } 

  }


  partial_C *= 1./this -> scaleFactor;
  partial_S *= 1./this -> scaleFactor;



}



void SBGATSphericalHarmo::PrintHeader(ostream& os, vtkIndent indent) {

}
void SBGATSphericalHarmo::PrintTrailer(ostream& os, vtkIndent indent) {

}


void SBGATSphericalHarmo::SaveToJson(std::string path) const{


  // The spherical harmonics are saved to a JSON file
  // The JSON fieds are:
  // - facets == number of facets
  // - vertices == number of vertices
  // - totalMass : {value, unit}
  // - density : {value, unit}
  // - reference_radius : {value, unit}
  // - normalized == true if the coefficients are normalized
  // - degree == degree of the spherical expansion
  // - Cnm_coefs - vector of coefficients triplets
  // - Snm_coefs - vector of coefficients triplets



  nlohmann::json spherical_harmo_json;

  spherical_harmo_json["facets"] = this -> n_facets;
  spherical_harmo_json["vertices"] = this -> n_vertices;
  spherical_harmo_json["totalMass"]["value"] = this -> totalMass;
  spherical_harmo_json["totalMass"]["unit"] = "kg";


  spherical_harmo_json["density"]["value"] = this -> density;

  spherical_harmo_json["referenceRadius"]["value"] = this -> referenceRadius;
  spherical_harmo_json["density"]["unit"] = "kg/m^3";

  if (this -> scaleFactor == 1){
    spherical_harmo_json["referenceRadius"]["unit"] = "m";

  }
  else{
    spherical_harmo_json["referenceRadius"]["unit"] = "km";

  }

  spherical_harmo_json["normalized"] = this -> normalized;
  spherical_harmo_json["degree"] = this -> degree;

  nlohmann::json Cnm_coefs;
  nlohmann::json Snm_coefs;


  for (unsigned int nn = 0; nn<=this -> degree; nn++){
    for (unsigned int mm = 0; mm<=nn; mm++){

      nlohmann::json coef_C = { {"n", nn}, {"m", mm}, {"value", this -> Cnm(nn,mm)} };
      nlohmann::json coef_S = { {"n", nn}, {"m", mm},{"value", this -> Snm(nn,mm)} };

      Cnm_coefs.push_back(coef_C);
      Snm_coefs.push_back(coef_S);
    }
  }

  spherical_harmo_json["Cnm_coefs"] = Cnm_coefs;
  spherical_harmo_json["Snm_coefs"] = Snm_coefs;


  std::ofstream o(path);
  o << std::setw(4) << spherical_harmo_json << std::endl;



}

void SBGATSphericalHarmo::LoadFromJson(std::string path){

  // The JSON container is created
  nlohmann::json spherical_harmo_json;

  // The file is loaded into the container
  std::ifstream i(path);  
  i >> spherical_harmo_json;

  // There should be a total of 9 objects in the container
  if (spherical_harmo_json.size() != 9){
    throw(std::runtime_error("Parsed JSON file should contain 9 objects, but SBGAT found " + std::to_string(spherical_harmo_json.size()) ));
  }

  // The fields are parsed and used to set the SBGATSphericalHarmo state
  this -> density = spherical_harmo_json.at("density").at("value");
  this -> referenceRadius = spherical_harmo_json.at("referenceRadius").at("value");
  this -> totalMass = spherical_harmo_json.at("totalMass").at("value");

  if (spherical_harmo_json.at("referenceRadius").at("unit") == "m" ){
    this -> scaleFactor = 1;
  }
  else{
    this -> scaleFactor = 1000;
  }

  this -> normalized = spherical_harmo_json.at("normalized");
  this -> degree = spherical_harmo_json.at("degree");

  this -> Cnm.clear();
  this -> Snm.clear();
  this -> Cnm = arma::zeros<arma::mat>(this -> degree + 1  , this -> degree + 1);
  this -> Snm = arma::zeros<arma::mat>(this -> degree + 1  , this -> degree + 1);


  nlohmann::json Cnm_coefs = spherical_harmo_json.at("Cnm_coefs");
  nlohmann::json Snm_coefs = spherical_harmo_json.at("Snm_coefs");


  for (nlohmann::json::iterator it = Cnm_coefs.begin(); it != Cnm_coefs.end(); ++it) {
    int n = it -> at("n");
    int m = it -> at("m");
    double value = it -> at("value");
    this -> Cnm(n,m) = value;
  }

  for (nlohmann::json::iterator it = Snm_coefs.begin(); it != Snm_coefs.end(); ++it) {
    int n = it -> at("n");
    int m = it -> at("m");
    double value = it -> at("value");
    this -> Snm(n,m) = value;
  }  

  this -> setFromJSON = true;

  // Will silence a warning thrown when Update() is called after the spherical harmonics 
  // were loaded from a JSON file since no vtkPolydata was effectively connected
  // to $this. Obviously $empty_polydata should not (and will not) be used 
  vtkSmartPointer<vtkPolyData> empty_polydata = vtkSmartPointer<vtkPolyData>::New();
  this -> SetInputData(empty_polydata);

}

//----------------------------------------------------------------------------
void SBGATSphericalHarmo::PrintSelf(std::ostream& os, vtkIndent indent){

  vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
  if (!input)
  {
    return;
  }

}
