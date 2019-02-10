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


#include "PickInteractorStyle.hpp"
#include <vtkRenderWindowInteractor.h>
#include <vtkCellPicker.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include <vtkSphereSource.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkLinearExtrusionFilter.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

#include "Mainwindow.hpp"

using namespace SBGAT_GUI;


vtkStandardNewMacro(PickInteractorStyle);


PickInteractorStyle::PickInteractorStyle(){
	this -> selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	this -> selectedActor = vtkSmartPointer<vtkActor>::New();
	this -> selectedActor->SetMapper(this -> selectedMapper);

}

void PickInteractorStyle::SetMainwindow(Mainwindow * mainwindow){
	this -> mainwindow = mainwindow;
}

void PickInteractorStyle::OnLeftButtonDown(){
	
	std::string name;
	// The picked shape corresponds to the selected one
	if (this -> mainwindow -> get_wrapped_shape_data().size() > 0){
		int selected_row_index = this -> mainwindow -> prop_table -> selectionModel() -> currentIndex().row();
		name = this ->mainwindow ->  prop_table -> item(selected_row_index, 0) -> text() .toStdString();
		this -> Data = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_polydata();
	}

    // Get the location of the click (in window coordinates)
	int * pos = this->GetInteractor()->GetEventPosition();

	vtkSmartPointer<vtkCellPicker> picker = vtkSmartPointer<vtkCellPicker>::New();
	picker -> SetTolerance(0.0005);

    // Pick from this location.
	picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());



	if(picker -> GetCellId() != -1){




		double bounds[6];
		this -> Data -> GetBounds(bounds);

		arma::vec mins = {bounds[0],bounds[2],bounds[4]};
		arma::vec maxs = {bounds[1],bounds[3],bounds[5]};

		this -> object_size = arma::abs(maxs - mins).max();

		if (this -> mainwindow -> get_selection_mode()){

			// Facet picking

			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
			triangle->GetPointIds()->SetId ( 0, 0 );
			triangle->GetPointIds()->SetId ( 1, 1 );
			triangle->GetPointIds()->SetId ( 2, 2);

			vtkCell * picked_cell = this -> Data ->GetCell(picker -> GetCellId());

			vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();

			double p0[3];
			double p1[3];
			double p2[3];
			picked_cell->GetPoints()->GetPoint(0, p0);
			picked_cell->GetPoints()->GetPoint(1, p1);
			picked_cell->GetPoints()->GetPoint(2, p2);

			vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
			points->InsertNextPoint ( p0[0], p0[1], p0[2] );
			points->InsertNextPoint ( p1[0], p1[1], p1[2] );
			points->InsertNextPoint ( p2[0], p2[1], p2[2] );


			vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
			triangles -> InsertNextCell ( triangle );

			polydata -> SetPoints ( points );
			polydata -> SetPolys ( triangles );

			double n[3];

			vtkTriangle::ComputeNormal (p0, p1, p2, n);

			vtkSmartPointer<vtkLinearExtrusionFilter> extrude = 
			vtkSmartPointer<vtkLinearExtrusionFilter>::New();
			extrude -> SetInputData( polydata);
			extrude -> SetExtrusionTypeToNormalExtrusion();
			extrude -> SetVector(n[0], n[1], n[2] );
			extrude -> SetScaleFactor (1);
			extrude -> CappingOn();

			this -> selectedMapper -> SetInputConnection(extrude -> GetOutputPort());

			this -> selectedActor -> GetProperty() -> VertexVisibilityOn();
			this -> selectedActor -> GetProperty() -> EdgeVisibilityOn();

			this -> selectedActor -> GetProperty() -> SetEdgeColor(1,0,0);
			this -> selectedActor -> GetProperty() -> SetVertexColor(1,0,0);
			this -> selectedActor -> GetProperty() -> SetColor(1,0,0);


			this -> Interactor -> GetRenderWindow() -> GetRenderers()->GetFirstRenderer()->AddActor(this -> selectedActor);

			this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("Picked facet " + std::to_string(picker->GetCellId())));

			vtkSmartPointer<vtkFloatArray> slopes = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_slopes();
			if (slopes != nullptr ){

				if(slopes -> GetNumberOfValues()> 0){

					auto slopes = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_slopes();
					auto inertial_potentials = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_inertial_potentials();
					auto body_fixed_potentials = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_body_fixed_potentials();
					auto inertial_acc_magnitudes = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_inertial_acc_magnitudes();
					auto body_fixed_acc_magnitudes = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_body_fixed_acc_magnitudes();

					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\tSlope: " + std::to_string(slopes->GetValue(picker->GetCellId())) + " (deg)"));
					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\tInertial potential: " + std::to_string(inertial_potentials->GetValue(picker->GetCellId())) + " (m^2/s^2)"));
					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\tBody-fixed potential: " + std::to_string(body_fixed_potentials->GetValue(picker->GetCellId())) + " (m^2/s^2)"));
					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\tInertial acceleration magnitude: " + std::to_string(inertial_acc_magnitudes->GetValue(picker->GetCellId())) + " (m/s^2)"));
					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\tBody-fixed acceleration magnitude: " + std::to_string(body_fixed_acc_magnitudes->GetValue(picker->GetCellId())) + " (m/s^2)"));
					this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("\n"));

				}
			}

		}
		else{

			// Points picking
			double * pick_coordinates = picker->GetPickPosition();
			this -> selected_point_id = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_tree() -> FindClosestPoint(pick_coordinates);

			double vertex_coordinates[3]; 
			vtkTriangle* triangle  = dynamic_cast<vtkTriangle*>(this -> Data ->GetCell(picker -> GetCellId()));

			double p0[3];
			double p1[3];
			double p2[3];
			triangle->GetPoints()->GetPoint(0, p0);
			triangle->GetPoints()->GetPoint(1, p1);
			triangle->GetPoints()->GetPoint(2, p2);
			double area = vtkTriangle::TriangleArea(p0, p1, p2);

			this -> Data -> GetPoint(this -> selected_point_id, vertex_coordinates); 

			this -> mainwindow -> log_console -> appendPlainText(QString::fromStdString("Picked vertex " + std::to_string(this -> selected_point_id)));
			
			// Create the geometry of a point (the coordinate)

			vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
			sphere -> SetCenter (vertex_coordinates[0], vertex_coordinates[1],  vertex_coordinates[2]);

			sphere -> SetPhiResolution (20);
			sphere -> SetThetaResolution (20);
			sphere -> SetRadius (0.25 * std::sqrt(area));


			sphere -> Update();

			this -> selectedMapper -> SetInputConnection(sphere -> GetOutputPort());
			this -> selectedActor -> GetProperty() -> SetColor(1,0,0);
			this -> selectedActor -> GetProperty() -> VertexVisibilityOff();
			this -> selectedActor -> GetProperty() -> EdgeVisibilityOff();



			this -> Interactor -> GetRenderWindow() -> GetRenderers() -> GetFirstRenderer() -> AddActor(this -> selectedActor);

			this -> GetInteractor() -> GetRenderWindow() -> Render();


		}

		
	}
	else{
		this -> Clear();
	}
	// Forward events
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void PickInteractorStyle::Clear(){
	this -> Interactor -> GetRenderWindow() -> GetRenderers() -> GetFirstRenderer() -> RemoveActor(this -> selectedActor);
	this -> Data = nullptr;
	this -> GetInteractor() -> GetRenderWindow() -> Render();
}		


int PickInteractorStyle::GetSelectionSize() const{

	if (this -> Data == nullptr){
		return 0;
	}
	return std::max(this -> selectedMapper -> GetInput() -> GetNumberOfPoints(),
		this -> selectedMapper -> GetInput() -> GetNumberOfCells());

}

void PickInteractorStyle::GetNormalAtSelectedPoint(double * normal) const{
	
	vtkFloatArray * normalDataFloat = vtkFloatArray::SafeDownCast(this -> Data ->GetPointData()->GetArray("Normals"));
	normalDataFloat -> GetTuple(this -> selected_point_id,normal);

}

void PickInteractorStyle::GetNormalAtPoint(vtkIdType id, double * normal) const{

	vtkFloatArray * normalDataFloat = vtkFloatArray::SafeDownCast(this -> Data ->GetPointData()->GetArray("Normals"));
	normalDataFloat -> GetTuple(id,normal);

}


void PickInteractorStyle::GetQueriedPoint(double * query){
	this -> Data -> GetPoint(this -> selected_point_id, query);

}
