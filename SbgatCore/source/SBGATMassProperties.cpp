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
  Module:    SBGATMassProperties.cxx

  Derived class from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "SBGATMassProperties.hpp"

#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>
vtkStandardNewMacro(SBGATMassProperties);

//----------------------------------------------------------------------------
// Constructs with initial 0 values.
SBGATMassProperties::SBGATMassProperties(){
  this->SurfaceArea = 0.0;
  this->MinCellArea = 0.0;
  this->MaxCellArea = 0.0;
  this->Volume  = 0.0;
  this->VolumeProjected = 0.0;
  this->VolumeX = 0.0;
  this->VolumeY = 0.0;
  this->VolumeZ = 0.0;
  this->Kx = 0.0;
  this->Ky = 0.0;
  this->Kz = 0.0;
  this->NormalizedShapeIndex = 0.0;

  this->SetNumberOfOutputPorts(0);
}

//----------------------------------------------------------------------------
// Destroy any allocated memory.
SBGATMassProperties::~SBGATMassProperties(){
}

//----------------------------------------------------------------------------
// Description:
// This method measures volume, surface area, normalized shape index, the surface area, the
//  normalized shape index, center of mass (assuming constant density) and inertia tensor
// (assuming constant density) of a triangle mesh.

int SBGATMassProperties::RequestData(
  vtkInformation* vtkNotUsed( request ),
  vtkInformationVector** inputVector,
  vtkInformationVector* vtkNotUsed( outputVector )){
  vtkInformation *inInfo =
  inputVector[0]->GetInformationObject(0);

  // call ExecuteData
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType cellId, numCells, numPts, numIds;
  double p[3];

  numCells = input->GetNumberOfCells();
  numPts = input->GetNumberOfPoints();
  if (numCells < 1 || numPts < 1)
  {
    vtkErrorMacro( << "No data to measure...!");
    return 1;
  }

  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  ptIds -> Allocate(VTK_CELL_SIZE);

  // Traverse all cells, obtaining node coordinates.
  //
  double    vol[3],kxyz[3];
  double    xp[3]; // to compute volumeproj
  double    munc[3],wxyz,wxy,wxz,wyz;
  double    area,surfacearea;
  double sum_surface[3];
  double    volumeproj;
  double    mincellarea, maxcellarea;
  double    a,b,c,s;
  double    x[3],y[3],z[3];
  double    i[3],j[3],k[3],u[3],absu[3],length;
  double    ii[3],jj[3],kk[3];
  double    xavg,yavg,zavg;
  double cx,cy,cz;
  double dv;
  double P_xx;
  double P_yy;
  double P_zz;
  double P_xy;
  double P_xz;
  double P_yz;
  double l; //normalizing length
  double average_surface;

  vtkIdType idx;

  // Initialize variables ...
  //
  surfacearea = 0.0;
  volumeproj = 0.0;
  mincellarea = VTK_DOUBLE_MAX; maxcellarea = 0.0;
  wxyz = 0; wxy = 0.0; wxz = 0.0; wyz = 0.0;
  cx = 0; cy = 0; cz = 0;
  P_xx = 0;   P_yy = 0;   P_zz = 0;
  P_xy = 0;   P_yz = 0;   P_xz = 0;

  sum_surface[0] = 0;
  sum_surface[1] = 0;
  sum_surface[2] = 0;
  average_surface = 0;


  input -> GetBounds(this -> bounds);
  arma::vec::fixed<3> bbox_min = {this -> bounds[0],this -> bounds[2],this -> bounds[4]};
  arma::vec::fixed<3> bbox_max = {this -> bounds[1],this -> bounds[3],this -> bounds[5]};
  l = arma::norm(bbox_max - bbox_min);

  for ( idx = 0; idx < 3 ; idx++ ){
    munc[idx] = 0.0;
    vol[idx]  = 0.0;
    kxyz[idx] = 0.0;
  }

  for (cellId=0; cellId < numCells; cellId++){
    if ( input->GetCellType(cellId) != VTK_TRIANGLE){
      vtkWarningMacro(<< "Input data type must be VTK_TRIANGLE not "<< input->GetCellType(cellId));
      continue;
    }
    input->GetCellPoints(cellId,ptIds);
    numIds = ptIds->GetNumberOfIds();
    assert(numIds == 3);

    // store current vertex (x,y,z) coordinates ...
    // Note that the coordinates are normalized!
    for (idx=0; idx < numIds; idx++){
      input->GetPoint(ptIds->GetId(idx), p);
      x[idx] = p[0] / l; y[idx] = p[1] / l; z[idx] = p[2] / l;
    }

    // get i j k vectors ...
    //
    i[0] = ( x[1] - x[0]); j[0] = (y[1] - y[0]); k[0] = (z[1] - z[0]);
    i[1] = ( x[2] - x[0]); j[1] = (y[2] - y[0]); k[1] = (z[2] - z[0]);
    i[2] = ( x[2] - x[1]); j[2] = (y[2] - y[1]); k[2] = (z[2] - z[1]);

    // cross product between two vectors, to determine normal vector
    //
    u[0] = ( j[0] * k[1] - k[0] * j[1]);
    u[1] = ( k[0] * i[1] - i[0] * k[1]);
    u[2] = ( i[0] * j[1] - j[0] * i[1]);

    // normalize normal vector to 1
    //
    length = sqrt( u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    if ( length != 0.0)
    {
      u[0] /= length;
      u[1] /= length;
      u[2] /= length;
    }
    else
    {
      u[0] = u[1] = u[2] = 0.0;
    }

    // determine max unit normal component...
    //
    absu[0] = fabs(u[0]); absu[1] = fabs(u[1]); absu[2] = fabs(u[2]);

    if (( absu[0] > absu[1]) && ( absu[0] > absu[2]) )
    {
      munc[0]++;
    }
    else if (( absu[1] > absu[0]) && ( absu[1] > absu[2]) )
    {
      munc[1]++;
    }
    else if (( absu[2] > absu[0]) && ( absu[2] > absu[1]) )
    {
      munc[2]++;
    }
    else if (( absu[0] == absu[1])&& ( absu[0] == absu[2]))
    {
      wxyz++;
    }
    else if (( absu[0] == absu[1])&& ( absu[0] > absu[2]) )
    {
      wxy++;
    }
    else if (( absu[0] == absu[2])&& ( absu[0] > absu[1]) )
    {
      wxz++;
    }
    else if (( absu[1] == absu[2])&& ( absu[0] < absu[2]) )
    {
      wyz++;
    }
    else
    {
      vtkErrorMacro( << "Unpredicted situation...!" );
      return 1;
    }

    // This is reduced to ...
    //
    ii[0] = i[0] * i[0]; ii[1] = i[1] * i[1]; ii[2] = i[2] * i[2];
    jj[0] = j[0] * j[0]; jj[1] = j[1] * j[1]; jj[2] = j[2] * j[2];
    kk[0] = k[0] * k[0]; kk[1] = k[1] * k[1]; kk[2] = k[2] * k[2];

    // area of a triangle...
    //
    a = sqrt(ii[1] + jj[1] + kk[1]);
    b = sqrt(ii[0] + jj[0] + kk[0]);
    c = sqrt(ii[2] + jj[2] + kk[2]);
    s = 0.5 * (a + b + c);
    area = sqrt( fabs(s*(s-a)*(s-b)*(s-c)));
    surfacearea += area;

    average_surface += area / numCells;



    if( area < mincellarea )
    {
      mincellarea = area;
    }
    if( area > maxcellarea )
    {
      maxcellarea = area;
    }

    // volume elements ...
    //
    zavg = (z[0] + z[1] + z[2]) / 3.0;
    yavg = (y[0] + y[1] + y[2]) / 3.0;
    xavg = (x[0] + x[1] + x[2]) / 3.0;

    vol[2] += (area * u[2] * zavg);
    vol[1] += (area * u[1] * yavg);
    vol[0] += (area * u[0] * xavg);

    // V  =  (z1+z2+z3)(x1y2-x2y1+x2y3-x3y2+x3y1-x1y3)/6
    // Volume under triangle is projected area of the triangle times
    // the average of the three z values
    vtkMath::Cross(x,y,xp);
    volumeproj += zavg * (xp[0]+xp[1]+xp[2]) / 2;

    // Center of mass
    // See "Inertia of Any Polyhedron" by Anthony R. Dobrovolskis, Icarus 124, 698–704 (1996) Article No. 0243
    dv = 1. / 6. * (x[1] * ( (  y[1] - y[0] ) * (  z[2] - z[0]   ) - (  z[1] - z[0]   ) * (  y[2] - y[0])) + 
      y[1] * ( (  z[1] - z[0] ) * (  x[2] - x[0]   ) -  (  x[1] - x[0]   ) * (  z[2] - z[0]   ) ) + 
      z[1] * ( (  x[1] - x[0] ) * (  y[2] - y[0]   ) -  (  y[1] - y[0]   ) * (  x[2] - x[0]   )));

    cx += dv * (x[0] + x[1] + x[2])/4;
    cy += dv * (y[0] + y[1] + y[2])/4;
    cz += dv * (z[0] + z[1] + z[2])/4;


    // Inertia tensor
    // See "Inertia of Any Polyhedron" by Anthony R. Dobrovolskis, Icarus 124, 698–704 (1996) Article No. 0243

    P_xx += dv / 20 * (2 * x[0] * x[0]
      + 2 * x[1] * x[1]
      + 2 * x[2] * x[2]
      + x[0] * x[1]
      + x[0] * x[1]
      + x[0] * x[2]
      + x[0] * x[2]
      + x[1] * x[2]
      + x[1] * x[2]);


    P_yy += dv / 20 * (2 * y[0] * y[0]
      + 2 * y[1] * y[1]
      + 2 * y[2] * y[2]
      + y[0] * y[1]
      + y[0] * y[1]
      + y[0] * y[2]
      + y[0] * y[2]
      + y[1] * y[2]
      + y[1] * y[2]);

    P_zz += dv / 20 * (2 * z[0] * z[0]
      + 2 * z[1] * z[1]
      + 2 * z[2] * z[2]
      + z[0] * z[1]
      + z[0] * z[1]
      + z[0] * z[2]
      + z[0] * z[2]
      + z[1] * z[2]
      + z[1] * z[2]);

    P_xy += dv / 20 * (2 * x[0] * y[0]
      + 2 * x[1] * y[1]
      + 2 * x[2] * y[2]
      + x[0] * y[1]
      + y[0] * x[1]
      + x[0] * y[2]
      + y[0] * x[2]
      + x[1] * y[2]
      + y[1] * x[2]);

    P_xz += dv / 20 * (2 * x[0] * z[0]
      + 2 * x[1] * z[1]
      + 2 * x[2] * z[2]
      + x[0] * z[1]
      + z[0] * x[1]
      + x[0] * z[2]
      + z[0] * x[2]
      + x[1] * z[2]
      + z[1] * x[2]);

    P_yz += dv / 20 * (2 * y[0] * z[0]
      + 2 * y[1] * z[1]
      + 2 * y[2] * z[2]
      + y[0] * z[1]
      + z[0] * y[1]
      + y[0] * z[2]
      + z[0] * y[2]
      + y[1] * z[2]
      + z[1] * y[2]);

    // Sum of oriented surface
    vtkMath::MultiplyScalar(u,area);
    vtkMath::Add(u,sum_surface,sum_surface);

  }

  // Surface Area ...
  //
  this->SurfaceArea = surfacearea* std::pow(l,2);
  this->MinCellArea = mincellarea* std::pow(l,2);
  this->MaxCellArea = maxcellarea* std::pow(l,2);

  // Weighting factors in Discrete Divergence theorem for volume calculation.
  //
  kxyz[0] = (munc[0] + (wxyz/3.0) + ((wxy+wxz)/2.0)) /numCells;
  kxyz[1] = (munc[1] + (wxyz/3.0) + ((wxy+wyz)/2.0)) /numCells;
  kxyz[2] = (munc[2] + (wxyz/3.0) + ((wxz+wyz)/2.0)) /numCells;
  this->VolumeX = vol[0] * std::pow(l,3) ;
  this->VolumeY = vol[1] * std::pow(l,3) ;
  this->VolumeZ = vol[2] * std::pow(l,3) ;
  this->Kx = kxyz[0];
  this->Ky = kxyz[1];
  this->Kz = kxyz[2];
  this->Volume =  (kxyz[0] * vol[0] + kxyz[1] * vol[1] + kxyz[2]  * vol[2])* std::pow(l,3);
  this->Volume =  fabs(this->Volume);
  this->VolumeProjected = volumeproj * std::pow(l,3);
  this->NormalizedShapeIndex =
  (sqrt(surfacearea)/std::cbrt(this->Volume))/2.199085233;


  // Center of mass
  arma::vec::fixed<3> com_ = {cx,cy,cz};

  this -> center_of_mass = l * com_ / this->Volume * std::pow(l,3);

  // Inertia tensor that was non-dimensionalized by l, not cbrt(V)
  arma::mat I = {
    {P_yy + P_zz, -P_xy, -P_xz},
    { -P_xy, P_xx + P_zz, -P_yz},
    { -P_xz, -P_yz, P_xx + P_yy}};

    double L = std::cbrt(this -> Volume) ;
    this -> inertia_tensor = std::pow(l / L,5) *  I - RBK::tilde(this -> center_of_mass / L) * RBK::tilde(this -> center_of_mass / L).t();

    // The principal axes are extracted
    arma::vec eig_val;
    arma::mat eig_vec;
    arma::eig_sym(eig_val,eig_vec,this -> inertia_tensor);


    // At this stage, eig_vec is (+/- v0,+/- v1 , +/- v2) where |v0| = |v1| = |v2| = 1
    // with their corresponding values sorted in ascending order
    // But it is not guaranteed that det(eig_vec) = +1
    // The 8 total configurations of the principal axes matrix M that exist are
    // M1: (v0,v1,v2) == eig_vec
    // M2: (v0,v1,-v2)
    // M3: (v0,-v1,v2)
    // M4: (v0,-v1,-v2)
    // M5: (-v0,v1,v2)
    // M6: (-v0,v1,-v2)
    // M7: (-v0,-v1,v2)
    // M8: (-v0,-v1,-v2)

    // Note that M1,M4,M6,M7 and M2,M3,M5,M8 have the same determinant respectively

    if (arma::det(eig_vec) < 0){
      eig_vec = - eig_vec;
    }

    if(eig_vec(0,0) > 0 && eig_vec(1,1) < 0){

      eig_vec.col(1) = - eig_vec.col(1);
      eig_vec.col(2) = - eig_vec.col(2);
    }

    else if(eig_vec(0,0) < 0 && eig_vec(2,2) < 0){

      eig_vec.col(0) = - eig_vec.col(0);
      eig_vec.col(2) = - eig_vec.col(2);
    }

    else if(eig_vec(0,0) < 0 && eig_vec(1,1) < 0){

      eig_vec.col(0) = - eig_vec.col(0);
      eig_vec.col(1) = - eig_vec.col(1);
    }

    this -> principal_axes = eig_vec;


    // Closeness of topology given sum of oriented surface
    if (vtkMath::Norm(sum_surface) / average_surface < 1e-6){
      this -> IsClosed = true;
    }
    else{
      this -> IsClosed = false;
    }

    return 1;
  }



  void SBGATMassProperties::PrintHeader(ostream& os, vtkIndent indent) {

  }
  void SBGATMassProperties::PrintTrailer(ostream& os, vtkIndent indent) {

  }


//----------------------------------------------------------------------------
  void SBGATMassProperties::PrintSelf(std::ostream& os, vtkIndent indent){

    vtkPolyData *input = vtkPolyData::SafeDownCast(this->GetInput(0));
    if (!input){
      return;
    }
    os << "\tVolumeX: " << this->GetVolumeX () << "\n";
    os << "\tVolumeY: " << this->GetVolumeY () << "\n";
    os << "\tVolumeZ: " << this->GetVolumeZ () << "\n";
    os << "\tKx: " << this->GetKx () << "\n";
    os << "\tKy: " << this->GetKy () << "\n";
    os << "\tKz: " << this->GetKz () << "\n";
    os << "\tVolume:  " << this->GetVolume  () << "\n";
  //os << indent << "Volume Projected:  " << this->GetVolumeProjected  () << "\n";
  //os << indent << "Volume Error:  " <<
  //  fabs(this->GetVolume() - this->GetVolumeProjected())   << "\n";
    os << "\tSurface Area: " << this->GetSurfaceArea () << "\n";
    os << "\tMin Cell Area: " << this->GetMinCellArea () << "\n";
    os << "\tMax Cell Area: " << this->GetMaxCellArea () << "\n";
    os << "\tNormalized Shape Index: "
    << this->GetNormalizedShapeIndex () << "\n";
  }
