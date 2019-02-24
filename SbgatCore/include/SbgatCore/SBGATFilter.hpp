/*=========================================================================

  Program:   Small Body Geophysical Analysis
  Module:    SBGATFilter.hpp

  Class derived from VTK's vtkPolyDataAlgorithm by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
@file SBGATFilter.hpp
@class  SBGATFilter
@author Benjamin Bercovici
@date February 2019

@brief A generic class serving a a shape container and providing a few useful methods
to access the shape's underlying facets, vertices and edges
Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
@copyright MIT License, Benjamin Bercovici and Jay McMahon, (c) Ken Martin, Will Schroeder, Bill Lorensen
*/

#ifndef SBGATFilter_hpp
#define SBGATFilter_hpp

#include <vtkFiltersCoreModule.h> // For export macro
#include <vtkPolyDataAlgorithm.h>
#include <armadillo>


class VTKFILTERSCORE_EXPORT SBGATFilter : public vtkPolyDataAlgorithm{
public:

  /**
   * Constructs with initial values of zero.
   */

  static SBGATFilter *New();

  vtkTypeMacro(SBGATFilter,vtkPolyDataAlgorithm);
  void PrintSelf(std::ostream& os, vtkIndent indent) override;
  void PrintHeader(std::ostream& os, vtkIndent indent) override;
  void PrintTrailer(std::ostream& os, vtkIndent indent) override;

  /**
  Sets the scale factor to 1, indicative that the polydata has its coordinates expressed in meters
  */
  void SetScaleMeters() { this -> scaleFactor = 1; this -> scaleFactorSet = true;}

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; this -> scaleFactorSet = true;}

  /**
  Returns the shape's scale factor
  @return scale factor
  */
  double GetScaleFactor() const {return this -> scaleFactor;}

  /**
  Sets polyhedron density
  @param density bulk density of polyhedron (kg/m^3)
  */
  void SetDensity(const double density){
    this -> density = density;
    this -> densitySet = true;
  }


  /**
  Get polyhedron density
  @return density bulk density of polyhedron (kg/m^3)
  */
  double GetDensity() const{
    return this -> density ;
  }

 
  /**
  Returns coordinates of the three vertices forming a facet
  @param[in] f facet index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  @param[out] r2 first vertex coordinates
  */
  void GetVerticesInFacet(const int & f,double * r0,double * r1, double * r2) const;

  /**
  Returns coordinates of point at center of facet
  @param f facet index
  @return coordinates of f-th facet center
  */
  arma::vec::fixed<3> GetFacetCenter(const int & f) const;

  /**
  Returns coordinates of the two vertices forming an edge
  @param[in] e edge index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  */
  void GetVerticesOnEdge(const int & e,double * r0,double * r1) const;


  /**
  Returns normal of facet f
  @param[in] f facet index
  @param[out] n normal at facet
  */
  void GetFacetNormal(const int & f, double * n) const{n[0] = this -> facet_normals[f][0];n[1] = this -> facet_normals[f][1];n[2] = this -> facet_normals[f][2];}


  /**
  Returns the indices of the two facets adjacent to the specified edge
  @param[in] e edge index
  @param[out] f0 index of first facet
  @param[out] f1 index of second facet
  */
  void GetIndicesOfAdjacentFacets(const int & e,int & f0, int & f1) const;


  /**
  Returns the indices of the vertices forming the prescribed edge
  @param[in] e edge index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  */
  void GetIndicesVerticesOnEdge(const int & e, int & v0,int & v1){v0 = this -> edges[e][0];v1 = this -> edges[e][1];}

  /**
  Returns the indices of the vertices forming the prescribed facet
  @param[in] f facet index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  @param[out] v2 third vertex index
  */
  void GetIndicesVerticesInFacet(const int & f, int & v0,int & v1,int & v2){v0 = this -> facets[f][0];v1 = this -> facets[f][1];v2 = this -> facets[f][2];}


  /**
  Returns the non-normalized facet normal at this facet
  @param f facet index
  @return non-normalized facet normal
  */
  arma::vec::fixed<3> GetNonNormalizedFacetNormal(const int & f) const;

   /**
  * Compute and return the bounding box (xmin,xmax,ymin,ymax,zmin,zmax) (m)
  */
  double * GetBoundingBox(){
    this -> Update(); return this -> bounds;
  }

   /**
  * Compute and return the bounding box (xmin,xmax,ymin,ymax,zmin,zmax) (m)
  */
  void GetBoundingBox(double * bounds_){
    this -> Update(); bounds_ = this -> bounds;
  }

    /**
  * Compute and return the bounding box (xmin,xmax,ymin,ymax,zmin,zmax) (m)
  */
  void GetBoundingBox(double & xmin,double & xmax,double & ymin,double & ymax,double & zmin,double & zmax){
    this -> Update(); 
    xmin = bounds[0];
    xmax = bounds[1];
    ymin = bounds[2];
    ymax = bounds[3];
    zmin = bounds[4];
    zmax = bounds[5];
  }


protected:
  SBGATFilter();
  ~SBGATFilter() override;

  
  /**
  Returns the connectivity table associated with vector Tf
  @param f facet index
  @return connectivity table
  */
  arma::sp_mat  PartialTfPartialC(const int & f) const;

 

  void Clear();

  int RequestData(vtkInformation* request,
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector) override;

  double ** facet_normals;
  double ** vertices;

  double scaleFactor = 1;
  double density;

  int ** edges;
  int ** facets;
  int ** edge_facets_ids;

  bool scaleFactorSet;
  bool densitySet;

  int N_facets;
  int N_edges;
  int N_vertices;

  double bounds[6];


private:
  SBGATFilter(const SBGATFilter&) = delete;
  void operator=(const SBGATFilter&) = delete;
};

#endif


