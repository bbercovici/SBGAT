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
#include "SBGATMassProperties.hpp"
#include "SBGATFilter.hpp"


#include <vtkObjectFactory.h>
#include <vtkCell.h>
#include <vtkDataObject.h>
#include <vtkIdList.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <RigidBodyKinematics.hpp>
#include <json.hpp>
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
  SBGATFilter::Clear();
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
  vtkInformation * r = nullptr;
  vtkInformationVector * o = nullptr;

  SBGATFilter::RequestData( r ,inputVector,o);
  
  // Traverse all cells, obtaining node coordinates.
  //
  double    vol[3],kxyz[3];
  double tempVolume;
  double    xp[3]; // to compute volumeproj
  double    munc[3],wxyz,wxy,wxz,wyz;
  double    area,surfacearea;
  double sum_surface[3];
  double    volumeproj;
  double    mincellarea, maxcellarea;
  double    a,b,c,s;
  double    i[3],j[3],k[3],u[3],absu[3],length;
  double    ii[3],jj[3],kk[3];
  double    xavg,yavg,zavg;
  double average_surface;
  arma::vec::fixed<3> tempCOM = {0,0,0};
  arma::mat::fixed<3,3> I_temp = arma::zeros<arma::mat>(3,3);

  // Initialize variables ...
  surfacearea = 0.0;
  volumeproj = 0.0;
  tempVolume = 0;
  mincellarea = VTK_DOUBLE_MAX; maxcellarea = 0.0;
  wxyz = 0; wxy = 0.0; wxz = 0.0; wyz = 0.0;
  
  sum_surface[0] = 0;
  sum_surface[1] = 0;
  sum_surface[2] = 0;
  average_surface = 0;


  for ( int idx = 0; idx < 3 ; idx++ ){
    munc[idx] = 0.0;
    vol[idx]  = 0.0;
    kxyz[idx] = 0.0;
  }

  for (int f = 0; f < this -> N_facets; f++){

    double    x[3],y[3],z[3];
    
    // Note that the coordinates are scaled to meters if need be!

    x[0] = this -> scaleFactor * this -> vertices[facets[f][0]][0] ; y[0] = this -> scaleFactor * this -> vertices[facets[f][0]][1] ; z[0] = this -> scaleFactor * this -> vertices[facets[f][0]][2] ;
    x[1] = this -> scaleFactor * this -> vertices[facets[f][1]][0] ; y[1] = this -> scaleFactor * this -> vertices[facets[f][1]][1] ; z[1] = this -> scaleFactor * this -> vertices[facets[f][1]][2] ;
    x[2] = this -> scaleFactor * this -> vertices[facets[f][2]][0] ; y[2] = this -> scaleFactor * this -> vertices[facets[f][2]][1] ; z[2] = this -> scaleFactor * this -> vertices[facets[f][2]][2] ;

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

    average_surface += area / this -> N_facets;

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
    double dv = this -> GetDeltaV(f) ;
    tempCOM += this -> GetDeltaCM(f) * dv;
    I_temp +=  this -> GetDeltaIOverDeltaV(f) * dv;
    tempVolume += dv;

    // Sum of oriented surface
    vtkMath::MultiplyScalar(u,area);
    vtkMath::Add(u,sum_surface,sum_surface);


  }

  // Surface Area ...
  this->SurfaceArea = surfacearea;
  this->MinCellArea = mincellarea;
  this->MaxCellArea = maxcellarea;

  // Weighting factors in Discrete Divergence theorem for volume calculation.
  //
  kxyz[0] = (munc[0] + (wxyz/3.0) + ((wxy+wxz)/2.0)) /this -> N_facets;
  kxyz[1] = (munc[1] + (wxyz/3.0) + ((wxy+wyz)/2.0)) /this -> N_facets;
  kxyz[2] = (munc[2] + (wxyz/3.0) + ((wxz+wyz)/2.0)) /this -> N_facets;

  this->VolumeX = vol[0] ;
  this->VolumeY = vol[1] ;
  this->VolumeZ = vol[2] ;
  this -> Volume = tempVolume;

  this->Kx = kxyz[0];
  this->Ky = kxyz[1];
  this->Kz = kxyz[2];
  this->VolumeProjected = volumeproj ;
  this->NormalizedShapeIndex =(sqrt(surfacearea)/std::cbrt(this->Volume))/2.199085233;

  // Center of mass
  this -> center_of_mass = tempCOM / this -> Volume ;


  this -> r_avg =  std::cbrt( 3./4. * this -> Volume / arma::datum::pi ) ;
  this -> inertia_tensor = I_temp / (this -> Volume * this -> r_avg * this -> r_avg);

    // The principal axes are extracted
  arma::vec eig_val;
  arma::mat eig_vec;

  arma::eig_sym(eig_val,eig_vec,this -> inertia_tensor);

  if (arma::det(eig_vec) < 0){
    eig_vec.col(2) *= -1;
  }

  arma::mat::fixed<3,3> PB_1 = eig_vec.t();
  arma::mat::fixed<3,3> PB_2,PB_3,PB_4;

  PB_2.row(0) = -PB_1.row(0);
  PB_2.row(1) = -PB_1.row(1);
  PB_2.row(2) = PB_1.row(2);

  PB_3.row(0) = PB_1.row(0);
  PB_3.row(1) = -PB_1.row(1);
  PB_3.row(2) = -PB_1.row(2);

  PB_4.row(0) = -PB_1.row(0);
  PB_4.row(1) = PB_1.row(1);
  PB_4.row(2) = -PB_1.row(2);

  arma::vec::fixed<4> mrp_norms;

  mrp_norms(0) = arma::norm(RBK::dcm_to_mrp(PB_1));
  mrp_norms(1) = arma::norm(RBK::dcm_to_mrp(PB_2));
  mrp_norms(2) = arma::norm(RBK::dcm_to_mrp(PB_3));
  mrp_norms(3) = arma::norm(RBK::dcm_to_mrp(PB_4));

  // The stable principal axes are extracted by finding the dcm that yields the 
  // smallest-norm MRP

  if (mrp_norms.index_min() == 0){
   this -> principal_axes = PB_1;
 }
 else if (mrp_norms.index_min() == 1){
   this -> principal_axes = PB_2;
 }
 else if (mrp_norms.index_min() == 2){
   this -> principal_axes = PB_3;
 }
 else{
   this -> principal_axes = PB_4;
 }



 this -> unit_density_inertia_tensor = this -> Volume * this -> r_avg * this -> r_avg * this -> inertia_tensor;
 this -> normalized_principal_moments = arma::eig_sym(this -> inertia_tensor);
 this -> unit_density_principal_moments = arma::eig_sym(this -> unit_density_inertia_tensor);


 this -> principal_dimensions = std::sqrt(5. / 2. / this -> Volume) * arma::sqrt(arma::vec::fixed<3>({
  unit_density_principal_moments(1) + unit_density_principal_moments(2) - unit_density_principal_moments(0),
  unit_density_principal_moments(0) + unit_density_principal_moments(2) - unit_density_principal_moments(1),
  unit_density_principal_moments(0) + unit_density_principal_moments(1) - unit_density_principal_moments(2)
}));



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



void SBGATMassProperties::ComputeAndSaveMassProperties(vtkSmartPointer<vtkPolyData> shape,std::string path){

  vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
  mass_properties -> SetInputData(shape);
  mass_properties -> SetScaleMeters();

  mass_properties -> Update();

  mass_properties -> SaveMassProperties(path);

}

void SBGATMassProperties::SaveMassProperties(std::string path) const {

  nlohmann::json mass_properties_json;
  std::string length_unit,surface_unit,volume_unit;

  length_unit = "m";
  surface_unit = "m^2";
  volume_unit = "m^3";



  nlohmann::json com_json = {
    {"value",{this -> center_of_mass(0), this -> center_of_mass(1), this -> center_of_mass(2)}},
    {"unit",length_unit}
  };

  nlohmann::json volume_json = {
    {"value", this -> Volume},
    {"unit",volume_unit}
  };


  nlohmann::json surface_area_json = {
    {"value", this -> SurfaceArea},
    {"unit",surface_unit}
  };

  nlohmann::json normalized_inertia_tensor_json = {
    {"value", {
      {this -> inertia_tensor(0,0),this -> inertia_tensor(0,1),this -> inertia_tensor(0,2)},
      {this -> inertia_tensor(1,0),this -> inertia_tensor(1,1),this -> inertia_tensor(1,2)},
      {this -> inertia_tensor(2,0),this -> inertia_tensor(2,1),this -> inertia_tensor(2,2)}}
    },
    {"unit","none"}
  };

  nlohmann::json  unit_density_inertia_tensor_json = {
    {"value", {
      {this -> unit_density_inertia_tensor(0,0),this -> unit_density_inertia_tensor(0,1),this -> unit_density_inertia_tensor(0,2)},
      {this -> unit_density_inertia_tensor(1,0),this -> unit_density_inertia_tensor(1,1),this -> unit_density_inertia_tensor(1,2)},
      {this -> unit_density_inertia_tensor(2,0),this -> unit_density_inertia_tensor(2,1),this -> unit_density_inertia_tensor(2,2)}}
    },
    {"unit","m^5"}
  };


  nlohmann::json principal_dimensions_json = {
    {"value",{this -> principal_dimensions(0),this -> principal_dimensions(1),this -> principal_dimensions(2)}},
    {"unit",length_unit}
  };

  nlohmann::json unit_density_principal_moments_json = {
    {"value",{this -> unit_density_principal_moments(0),this -> unit_density_principal_moments(1),this -> unit_density_principal_moments(2)}},
    {"unit","m^5"}
  };



  nlohmann::json normalized_principal_moments_json = {
    {"value",{this -> normalized_principal_moments(0),this -> normalized_principal_moments(1),this -> normalized_principal_moments(2)}},
    {"unit","none"}
  };


  nlohmann::json projected_volume_components_json = {
    {"value",{this->VolumeX,this->VolumeY,this->VolumeZ}},
    {"unit",volume_unit}
  };

  nlohmann::json projection_coefs_json = {
    {"value",{this->Kx,this->Ky,this->Kz}},
    {"unit","none"}
  };

  nlohmann::json volume_projected_json = {
    {"value",this->VolumeProjected},
    {"unit",volume_unit}
  };

  nlohmann::json normalized_shape_index_json = {
    {"value",this->NormalizedShapeIndex},
    {"unit","none"}
  };

  nlohmann::json is_open_json = {
    {"value",this -> IsClosed ? "True" : "False"}
  };

  nlohmann::json average_radius_json = {
    {"value",this -> r_avg},
    {"unit",length_unit}
  };

  mass_properties_json["VOLUME"] = volume_json;
  mass_properties_json["SURFACE_AREA"] = surface_area_json;
  mass_properties_json["COM"] = com_json;
  mass_properties_json["VOLUME_PROJECTED"] = volume_projected_json;
  mass_properties_json["PROJECTED_VOLUME_VX_VY_VZ"] = projected_volume_components_json;
  mass_properties_json["PROJECTED_VOLUME_COEFS_KX_KY_KZ"] = projection_coefs_json;
  mass_properties_json["NORMALIZED_SHAPE_INDEX_JSON"] = normalized_shape_index_json;
  mass_properties_json["IS_CLOSED"] = is_open_json;
  mass_properties_json["NORMALIZED_INERTIA_TENSOR"] = normalized_inertia_tensor_json;
  mass_properties_json["NORMALIZED_PRINCIPAL_MOMENTS"] = normalized_principal_moments_json;
  mass_properties_json["UNIT_DENSITY_INERTIA_TENSOR"] = unit_density_inertia_tensor_json;
  mass_properties_json["UNIT_DENSITY_PRINCIPAL_MOMENTS"] = unit_density_principal_moments_json;
  mass_properties_json["PRINCIPAL_DIMENSIONS"] = principal_dimensions_json;

  mass_properties_json["AVERAGE_RADIUS"] = average_radius_json;


  std::ofstream o(path);
  o << std::setw(4) << mass_properties_json << std::endl;


}

double SBGATMassProperties::GetDeltaV(const int & f) const {

  return 1./6 * vtkMath::Determinant3x3(
    this -> vertices[this -> facets[f][0]],
    this -> vertices[this -> facets[f][1]],
    this -> vertices[this -> facets[f][2]]) * std::pow(this -> scaleFactor,3);

}

arma::vec::fixed<3> SBGATMassProperties::GetDeltaCM(const int & f) const {

  return (this -> scaleFactor/4 * arma::vec::fixed<3> ({
    this -> vertices[this -> facets[f][0]][0] + this -> vertices[this -> facets[f][1]][0] + this -> vertices[this -> facets[f][2]][0],
    this -> vertices[this -> facets[f][0]][1] + this -> vertices[this -> facets[f][1]][1] + this -> vertices[this -> facets[f][2]][1],
    this -> vertices[this -> facets[f][0]][2] + this -> vertices[this -> facets[f][1]][2] + this -> vertices[this -> facets[f][2]][2]})
  );

}


arma::mat::fixed<3,3> SBGATMassProperties::GetDeltaIOverDeltaV(const int & f) const {

  double r0[3];
  double r1[3];
  double r2[3];

  this -> GetVerticesInFacet(f,r0,r1,r2);
  arma::vec::fixed<3> C0 = {r0[0],r0[1],r0[2]};
  arma::vec::fixed<3> C1 = {r1[0],r1[1],r1[2]};
  arma::vec::fixed<3> C2 = {r2[0],r2[1],r2[2]};


  arma::mat::fixed<3,3> I = 6./5 * (
    1./2 * RBK::tilde(C0) * RBK::tilde(C0)
    + 1./6 * (RBK::tilde(C0) * RBK::tilde(C1 - C0) + RBK::tilde(C1 - C0) * RBK::tilde(C0))
    + 1./6 * (RBK::tilde(C0) * RBK::tilde(C2 - C0) + RBK::tilde(C2 - C0) * RBK::tilde(C0))
    + 1./24 * (RBK::tilde(C1 - C0) * RBK::tilde(C2 - C0) + RBK::tilde(C2 - C0) * RBK::tilde(C1 - C0))
    + 1./12 * RBK::tilde(C1 - C0) * RBK::tilde(C1 - C0)
    + 1./12 * RBK::tilde(C2 - C0) * RBK::tilde(C2 - C0)
    );
  return - I * std::pow(this -> scaleFactor,2);

}



