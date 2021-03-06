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
  void SetScaleMeters() { this -> scaleFactor = 1;this -> is_in_meters = true;}

  /**
  Sets the scale factor to 1000, indicative that the polydata has its coordinates expressed in kilometers
  */
  void SetScaleKiloMeters() { this -> scaleFactor = 1000; this -> is_in_meters = false;}

  /**
  Return the shape's scale factor
  @return scale factor
  */
  double GetScaleFactor() const {return this -> scaleFactor;}

  /**
  Sets polyhedron density
  @param density bulk density of polyhedron (kg/m^3)
  */
  void SetDensity(const double density){
    this -> density = density;
  }


  /**
  Get polyhedron density
  @return density bulk density of polyhedron (kg/m^3)
  */
  double GetDensity() const{
    return this -> density ;
  }

 

  /**
  Return coordinates of the three vertices forming a facet
  @param[in] f facet index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  @param[out] r2 first vertex coordinates
  */
  void GetVerticesInFacet(const int & f,double * r0,double * r1, double * r2) const;

  /**
  Return coordinates of point at center of facet
  @param f facet index
  @return coordinates of f-th facet center
  */
  arma::vec::fixed<3> GetFacetCenter(const int & f) const;

  /**
  Return coordinates of the two vertices forming an edge
  @param[in] e edge index
  @param[out] r0 first vertex coordinates
  @param[out] r1 first vertex coordinates
  */
  void GetVerticesOnEdge(const int & e,double * r0,double * r1) const;


  /**
  Return normal of facet f
  @param[in] f facet index
  @param[out] n normal at facet
  */
  void GetFacetNormal(const int & f, double * n) const{n[0] = this -> facet_normals[f][0];n[1] = this -> facet_normals[f][1];n[2] = this -> facet_normals[f][2];}


  /**
  Return the indices of the two facets adjacent to the specified edge
  @param[in] e edge index
  @param[out] f0 index of first facet
  @param[out] f1 index of second facet
  */
  void GetIndicesOfAdjacentFacets(const int & e,int & f0, int & f1) const;


  /**
  Return the indices of the vertices forming the prescribed edge
  @param[in] e edge index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  */
  void GetIndicesVerticesOnEdge(const int & e, int & v0,int & v1) const {v0 = this -> edges[e][0];v1 = this -> edges[e][1];}

  /**
  Return the indices of the vertices forming the prescribed facet
  @param[in] f facet index
  @param[out] v0 first vertex index
  @param[out] v1 second vertex index
  @param[out] v2 third vertex index
  */
  void GetIndicesVerticesInFacet(const int & f, int & v0,int & v1,int & v2) const {v0 = this -> facets[f][0];v1 = this -> facets[f][1];v2 = this -> facets[f][2];}


  /**
  Return the non-normalized facet normal at this facet
  @param f facet index
  @return non-normalized facet normal
  */
  arma::vec::fixed<3> GetNonNormalizedFacetNormal(const int & f) const;

  

  /**
  Return the vector from the field point to the first vertex on the designated edge
  @param pos field point
  @param e edge index
  @return the vector from the field point to the first vertex on the designated edge
  */
  arma::vec::fixed<3> GetRe(const arma::vec::fixed<3> & pos,const int & e) const;


  /**
  Return the vector from the field point to the first vertex on the designated edge
  @param pos field point
  @param e edge index
  @return the vector from the field point to the first vertex on the designated edge
  */
  arma::vec::fixed<3> GetRe(const double * pos,const int & e) const;



  /**
  Return the vector from the field point to the first vertex on the designated facet
  @param pos field point
  @param f facet index
  @return the vector from the field point to the first vertex on the designated facet
  */
  arma::vec::fixed<3> GetRf(const arma::vec::fixed<3> & pos,const int & f) const;


  /**
  Return the vector from the field point to the first vertex on the designated facet
  @param pos field point
  @param f facet index
  @return the vector from the field point to the first vertex on the designated facet
  */
  arma::vec::fixed<3> GetRf(const double * pos,const int & f) const;


  /**
  * Return the bounding box (xmin,xmax,ymin,ymax,zmin,zmax) (m)
  @param[in] xmin minimum x-bound of the shape (m)
  @param[in] xmax maximum x-bound of the shape (m)
  @param[in] ymin minimum y-bound of the shape (m)
  @param[in] ymax maximum y-bound of the shape (m)
  @param[in] zmin minimum z-bound of the shape (m)
  @param[in] zmax maximum z-bound of the shape (m)
  */
  void GetBoundingBox(double & xmin,double & xmax,double & ymin,double & ymax,double & zmin,double & zmax) const {
    xmin = bounds[0];
    xmax = bounds[1];
    ymin = bounds[2];
    ymax = bounds[3];
    zmin = bounds[4];
    zmax = bounds[5];
  }


/**
  Return length of considered edge in meters
  @param e edge index
  @return edge length(m )
*/  
  double GetEdgeLength(const int & e) const;


  /**
  Return number of facets in shape
  @return number of facets
  */
  int GetN_facets() const{return this -> N_facets;}
  /**
  Return number of edges in shape
  @return number of edges
  */
  int GetN_edges() const{return this -> N_edges;}

  /**
  Return number of vertices in shape
  @return number of vertices
  */
  int GetN_vertices() const{return this -> N_vertices;}

  
  /**
  Return true if the input shape has its coordinates expressed in meters, false otherwise
  */
  bool GetIsInMeters() const{ return this -> is_in_meters;}


  /**
  Return true if the input shape has its coordinates expressed in kilometers, false otherwise
  */
  bool GetIsInKiloMeters() const{ return !this -> is_in_meters;}


protected:
  SBGATFilter();
  ~SBGATFilter() override;


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

  
  bool is_in_meters = true;

  int N_facets;
  int N_edges;
  int N_vertices;

  double bounds[6];


private:
  SBGATFilter(const SBGATFilter&) = delete;
  void operator=(const SBGATFilter&) = delete;
};

#endif


