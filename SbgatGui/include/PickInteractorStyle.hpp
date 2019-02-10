
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
  Module:    PickInteractorStyle.hpp

  Derived class from VTK examples by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef PICKINTERACTOR_STYLE_HEADER
#define PICKINTERACTOR_STYLE_HEADER

#include <vtkPolyData.h>
#include <vtkDataSetMapper.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleTrackballCamera.h>

namespace SBGAT_GUI {
  class Mainwindow;

// Catch mouse events
  class PickInteractorStyle : public vtkInteractorStyleTrackballCamera {
  public:
    static PickInteractorStyle* New();

    vtkTypeMacro(PickInteractorStyle,vtkInteractorStyleTrackballCamera);

    PickInteractorStyle();
    void SetMainwindow(Mainwindow * mainwindow);

    /**
    Returns the number of plotted entities in the current selection. Does not necessarily return 1 if one
    facet is selected, since the facet actor is constituted by a number of edges/points
    */
    int GetSelectionSize() const;

    void Clear();

    virtual void OnLeftButtonDown();

    /**
  Returns world coordinates of query point 
    */
    void GetQueriedPoint(double * query);

    /**
    Get coordinates of normal at queried point
    */
    void GetNormalAtSelectedPoint(double * normal) const;

     /**
    Get coordinates of normal at specified point
    */
    void GetNormalAtPoint(vtkIdType id, double * normal) const;

    /**
    Returns selected actor
    */
    vtkSmartPointer<vtkActor> GetSelectedActor(){return this -> selectedActor;}

    /**
  Returns index of query point 
    */
    vtkIdType GetQueriedPointId() const {return this -> selected_point_id;}


    double GetSelectedBodySize() const{return this -> object_size;}


    vtkSmartPointer<vtkPolyData> Data = nullptr;
    double object_size;
    vtkIdType selected_point_id;
    vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    vtkSmartPointer<vtkActor> selectedActor;
    Mainwindow * mainwindow;

  };
}

#endif