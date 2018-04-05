/** MIT License

Copyright (c) 2018 Benjamin Bercovici

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


#include "InteractorStyle.hpp"

// shortcut to interactor modes
#define INTERACTOR_IS_ORIENT 0
#define INTERACTOR_IS_SELECT 1


InteractorStyle::InteractorStyle() {
	this -> CurrentMode = INTERACTOR_IS_ORIENT;
	this -> mainwindow = nullptr;
}

void InteractorStyle::OnLeftButtonUp() {
	// Forward events
	vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

	if (this -> mainwindow != nullptr) {
		// If the interactor is in selection mode
		// AND if the selection widget is not already open

		if (this -> CurrentMode == INTERACTOR_IS_SELECT
		         &&  dynamic_cast<ModifyAreaWidget *>(this -> mainwindow -> right_dockwidget -> widget()) != nullptr) {
			if (dynamic_cast<ModifyAreaWidget *>(this -> mainwindow -> right_dockwidget -> widget()) -> is_set() == false) {

				this -> grab_area();

			}

		}

	}

}

vtkStandardNewMacro(InteractorStyle);

void InteractorStyle::grab_area() {
	InheritedPicker * picker = static_cast< InheritedPicker * >(this -> Interactor -> GetPicker());

	vtkSmartPointer<vtkHardwareSelector> selector =
	    vtkSmartPointer<vtkHardwareSelector>::New();
	selector -> SetRenderer(this -> mainwindow -> get_renderer());

	selector -> SetArea(&picker -> get_dimensions_uint()[0]);
	selector -> SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
	vtkSmartPointer<vtkSelection> selection = selector -> Select();

	// If the selector grabbed something
	if (selection -> GetNumberOfNodes() > 0) {

		vtkSmartPointer<vtkExtractSelection> extractSelection =
		    vtkSmartPointer<vtkExtractSelection>::New();

		extractSelection -> SetInputData(0, this -> mainwindow -> get_asteroid() -> get_polydata());
		extractSelection -> SetInputData(1, selection);
		extractSelection -> Update();

		vtkSmartPointer<vtkGeometryFilter> selected_geometry =
		    vtkSmartPointer<vtkGeometryFilter>::New();

		selected_geometry -> SetInputConnection(extractSelection -> GetOutputPort());
		selected_geometry -> Update();

		vtkSmartPointer<vtkPolyData> selected_polydata = selected_geometry -> GetOutput();



		// The selection is inverted to find the cells that were not selected
		vtkSmartPointer<vtkSelectionNode> selectionNode = selection -> GetNode(0);
		selectionNode -> SetFieldType(vtkSelectionNode::CELL);
		selectionNode -> SetContentType(vtkSelectionNode::INDICES);
		selectionNode -> SetSelectionList(selected_polydata -> GetCellData() -> GetArray("ids"));
		selectionNode -> GetProperties() -> Set(vtkSelectionNode::INVERSE(), 1);

		extractSelection -> Update();

		vtkSmartPointer<vtkGeometryFilter> unselected_geometry =
		    vtkSmartPointer<vtkGeometryFilter>::New();

		unselected_geometry -> SetInputConnection(extractSelection -> GetOutputPort());
		unselected_geometry -> Update();


		vtkSmartPointer<vtkPolyData> unselected_polydata = unselected_geometry -> GetOutput();

		// The global indices of the unselected points are
		// stored in a vtkIdList
		vtkSmartPointer<vtkDataArray> unselected_ids = unselected_polydata -> GetPointData() -> GetArray("ids");

		vtkSmartPointer<vtkIdList> unselected_vertex_ids_list = vtkSmartPointer<vtkIdList>::New();
		for (unsigned int i = 0; i < unselected_polydata -> GetNumberOfPoints(); ++i) {
			unselected_vertex_ids_list -> InsertNextId(*(unselected_ids -> GetTuple(i)));
		}


		// The same is done with the indices of the selected points
		vtkSmartPointer<vtkDataArray> selected_ids = selected_polydata -> GetPointData() -> GetArray("ids");
		vtkSmartPointer<vtkIdList> selected_vertex_ids_list = vtkSmartPointer<vtkIdList>::New();
		for (unsigned int i = 0; i < selected_polydata -> GetNumberOfPoints(); ++i) {
			selected_vertex_ids_list -> InsertNextId(*(selected_ids -> GetTuple(i)));
		}

		// A third list is created to store the boundary vertices.
		vtkSmartPointer<vtkIdList> boundary_vertex_ids_list = vtkSmartPointer<vtkIdList>::New();

		// This list is constructed by browsing through the shortest of the two lists
		// Vertices found in both the selected and unselectee polydata are on the boundary
		if (selected_polydata -> GetNumberOfPoints() < unselected_polydata -> GetNumberOfPoints()) {
			for (unsigned int i = 0; i < selected_polydata -> GetNumberOfPoints(); ++i) {

				if (unselected_vertex_ids_list -> IsId (selected_vertex_ids_list
				                                        -> GetId(i) ) != -1) {
					boundary_vertex_ids_list -> InsertNextId(selected_vertex_ids_list
					        -> GetId(i) );

				}
			}

		}
		else {

			for (unsigned int i = 0; i < unselected_polydata -> GetNumberOfPoints(); ++i) {

				if (selected_vertex_ids_list -> IsId (unselected_vertex_ids_list
				                                      -> GetId(i) ) != -1) {
					boundary_vertex_ids_list -> InsertNextId(unselected_vertex_ids_list
					        -> GetId(i) );

				}
			}

		}

		vtkSmartPointer<vtkPolyDataMapper> sel_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		sel_mapper -> SetInputData(selected_polydata);
		sel_mapper -> ScalarVisibilityOff();

		vtkSmartPointer<vtkActor> selected_actor = vtkSmartPointer<vtkActor>::New();
		selected_actor -> SetMapper(sel_mapper);
		selected_actor -> GetProperty() -> SetColor(1, 0, 0); //

		vtkSmartPointer<vtkPolyDataMapper> unsel_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
		unsel_mapper -> SetInputData(unselected_polydata);
		unsel_mapper -> ScalarVisibilityOff();

		vtkSmartPointer<vtkActor> unselected_actor = vtkSmartPointer<vtkActor>::New();
		unselected_actor -> SetMapper(unsel_mapper);

		this -> mainwindow -> get_renderer() -> AddActor(selected_actor);
		this -> mainwindow -> get_renderer() -> AddActor(unselected_actor);
		this -> mainwindow -> get_renderer() -> GetRenderWindow() -> Render();

		if (dynamic_cast<ModifyAreaWidget *>(this -> mainwindow -> right_dockwidget -> widget()) != nullptr) {


			dynamic_cast<ModifyAreaWidget *>(this -> mainwindow -> right_dockwidget -> widget()) -> set_data(
			    selected_polydata,
			    unselected_polydata,
			    selected_actor,
			    unselected_actor,
			    boundary_vertex_ids_list);
		}

		else {

			ModifyAreaWidget * modify_area_widget = new ModifyAreaWidget(
			    this -> mainwindow,
			    selected_polydata,
			    unselected_polydata,
			    selected_actor,
			    unselected_actor,
			    boundary_vertex_ids_list);

			this -> mainwindow -> right_dockwidget -> setWidget(modify_area_widget);
			this -> mainwindow -> right_dockwidget -> show();

			// the actors owned by the mainwindow are hidden. They will be
			// shown again once pc_editing widget is closed.
			// hiding them after highlighting the cells
			// prevents flickering
			this -> mainwindow -> set_actors_visibility(false);



		}



		// Back to orient mode
		this -> CurrentMode = INTERACTOR_IS_ORIENT;
	}
}



void InteractorStyle::set_mainwindow(Mainwindow * mainwindow) {
	this -> mainwindow = mainwindow;

}

void InteractorStyle::set_current_mode(const int mode) {
	this -> CurrentMode = mode;
}


Mainwindow * InteractorStyle::get_mainwindow() {
	return this -> mainwindow;
}


vtkSmartPointer<vtkPolyData> InteractorStyle::get_selected_points_polydata() {
	return this -> selected_points_polydata;
}