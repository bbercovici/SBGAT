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
		        &&  this -> mainwindow -> lateral_dockwidget -> widget() == nullptr) {

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

				// I'm not happy with the output of vtkFeatureEdges. For now,
				// it is possible to detect the edge by looking for vertices present in
				// both sets and constructing a boundary list this way

				// The boundary edges are found
				// vtkSmartPointer<vtkFeatureEdges> featureEdges =
				//     vtkSmartPointer<vtkFeatureEdges>::New();
				// featureEdges -> SetInputData(selected_polydata);
				// featureEdges -> BoundaryEdgesOn();
				// featureEdges -> FeatureEdgesOff();
				// featureEdges -> ManifoldEdgesOff();
				// featureEdges -> NonManifoldEdgesOff();
				// featureEdges -> Update();


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

				ModifyAreaWidget * modify_area_widget = new ModifyAreaWidget(
				    this -> mainwindow,
				    selected_polydata,
				    unselected_polydata,
				    selected_actor,
				    unselected_actor,
				    boundary_vertex_ids_list);

				this -> mainwindow ->  lateral_dockwidget -> setWidget(modify_area_widget);
				this -> mainwindow ->  lateral_dockwidget -> show();

				// the actors owned by the mainwindow are hidden. They will be
				// shown again once pc_editing widget is closed.
				// hiding them after highlighting the cells
				// prevents flickering
				this -> mainwindow -> set_actors_visibility(false);

				// Back to orient mode
				this -> CurrentMode = INTERACTOR_IS_ORIENT;

			}




		}

	}

}

vtkStandardNewMacro(InteractorStyle);



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