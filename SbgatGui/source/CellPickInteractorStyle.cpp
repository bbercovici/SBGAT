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
  Module:    CellPickInteractorStyle.hpp

  Derived class from VTK examples by Benjamin Bercovici  

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


#include "CellPickInteractorStyle.hpp"
#include <vtkRenderWindowInteractor.h>
#include <vtkCellPicker.h>
#include <vtkSelectionNode.h>
#include <vtkSelection.h>
#include <vtkExtractSelection.h>
#include <vtkUnstructuredGrid.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRendererCollection.h>
#include "Mainwindow.hpp"

using namespace SBGAT_GUI;


vtkStandardNewMacro(CellPickInteractorStyle);


CellPickInteractorStyle::CellPickInteractorStyle(){
	selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	selectedActor = vtkSmartPointer<vtkActor>::New();
}

void CellPickInteractorStyle::SetMainwindow(Mainwindow * mainwindow){
	this -> mainwindow = mainwindow;
}

void CellPickInteractorStyle::OnLeftButtonDown(){
	
	std::string name;
	// The picked shape corresponds to the selected one
	if (this -> mainwindow -> get_wrapped_shape_data().size() > 0){
		int selected_row_index = this -> mainwindow -> prop_table -> selectionModel() -> currentIndex().row();
		name = this ->mainwindow ->  prop_table -> item(selected_row_index, 0) -> text() .toStdString();
		this -> Data = this -> mainwindow -> get_wrapped_shape_data()[name] -> get_polydata();
	}


    // Get the location of the click (in window coordinates)
	int* pos = this->GetInteractor()->GetEventPosition();

	vtkSmartPointer<vtkCellPicker> picker =
	vtkSmartPointer<vtkCellPicker>::New();
	picker->SetTolerance(0.0005);

      // Pick from this location.
	picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());

	if(picker->GetCellId() != -1){

		
		vtkSmartPointer<vtkIdTypeArray> ids =
		vtkSmartPointer<vtkIdTypeArray>::New();
		ids->SetNumberOfComponents(1);
		ids->InsertNextValue(picker->GetCellId());

		vtkSmartPointer<vtkSelectionNode> selectionNode =
		vtkSmartPointer<vtkSelectionNode>::New();
		selectionNode->SetFieldType(vtkSelectionNode::CELL);
		selectionNode->SetContentType(vtkSelectionNode::INDICES);
		selectionNode->SetSelectionList(ids);

		vtkSmartPointer<vtkSelection> selection =
		vtkSmartPointer<vtkSelection>::New();
		selection->AddNode(selectionNode);

		vtkSmartPointer<vtkExtractSelection> extractSelection =
		vtkSmartPointer<vtkExtractSelection>::New();

		extractSelection->SetInputData(0, this->Data);
		extractSelection->SetInputData(1, selection);
		extractSelection->Update();

        // In selection
		vtkSmartPointer<vtkUnstructuredGrid> selected =
		vtkSmartPointer<vtkUnstructuredGrid>::New();
		selected->ShallowCopy(extractSelection->GetOutput());

		selectedMapper->SetInputData(selected);

		this -> selectedActor->SetMapper(selectedMapper);
		this -> selectedActor->GetProperty() -> EdgeVisibilityOn();
		this -> selectedActor->GetProperty() -> SetEdgeColor(1,0,0);
		this -> selectedActor->GetProperty() -> SetLineWidth(10);

		this -> Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(this -> selectedActor);

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
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this -> selectedActor);
	}
	// Forward events
	vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
}

void CellPickInteractorStyle::Clear(){
		this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->RemoveActor(this -> selectedActor);
}

