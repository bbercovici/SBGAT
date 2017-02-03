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

				// The boundary edges are found
				vtkSmartPointer<vtkFeatureEdges> featureEdges =
				    vtkSmartPointer<vtkFeatureEdges>::New();
				featureEdges -> SetInputData(selected_polydata);
				featureEdges -> BoundaryEdgesOn();
				featureEdges -> FeatureEdgesOff();
				featureEdges -> ManifoldEdgesOff();
				featureEdges -> NonManifoldEdgesOff();
				featureEdges -> Update();

				vtkSmartPointer<vtkDataArray> ids = featureEdges -> GetOutput() -> GetPointData() -> GetArray("ids");

				// The vertices located on the boundary of the polydata are located
				vtkSmartPointer<vtkIdList> boundary_vertex_ids_list = vtkSmartPointer<vtkIdList>::New();
				for (unsigned int i = 0; i < featureEdges -> GetOutput() -> GetNumberOfPoints(); ++i) {
					boundary_vertex_ids_list -> InsertNextId(*(ids -> GetTuple(i)));
					std::cout << *(ids -> GetTuple(i)) << std::endl;
				}
				std::cout << std::endl;

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


			// vtkSmartPointer<vtkSelectVisiblePoints> select_visible_points =
			//     vtkSmartPointer<vtkSelectVisiblePoints>::New();

			// vtkSmartPointer<vtkHardwareSelector> selector = vtkSmartPointer<vtkHardwareSelector>::New();
			// selector -> SetRenderer(this -> mainwindow -> get_renderer());

			// selector -> SetArea(&picker -> get_dimensions_uint()[0]);
			// selector -> SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
			// vtkSmartPointer<vtkSelection> selection = selector -> Select();


			// 	vtkSmartPointer<vtkExtractSelectedPolyDataIds> extractSelection =
			// 	    vtkSmartPointer<vtkExtractSelectedPolyDataIds>::New();

			// 	extractSelection -> SetInputData(0, this -> mainwindow -> get_asteroid() -> get_polydata());
			// 	extractSelection -> SetInputData(1, selection);
			// 	extractSelection -> Update();

			// 	vtkSmartPointer<vtkPolyData> selected_polydata_test = vtkPolyData::SafeDownCast(extractSelection -> GetOutput());
			// 	std::cout << selected_polydata_test -> GetNumberOfCells() << std::endl;
			// 	selection -> GetNode(0) -> Print(std::cout);
			// 	selection -> GetNode(0) -> GetProperties() -> Set(vtkSelectionNode::INVERSE(), 1); //invert the selection
			// 	selection -> GetNode(0) -> Print(std::cout);

			// 	extractSelection -> Update();

			// 	vtkSmartPointer<vtkPolyData> unselected_polydata_test = vtkPolyData::SafeDownCast(extractSelection -> GetOutput());
			// 	std::cout << unselected_polydata_test -> GetNumberOfCells() << std::endl;

			// 	// The filter renderer is set to the current renderer
			// 	select_visible_points -> SetRenderer(this -> mainwindow -> get_renderer());
			// 	select_visible_points -> SelectionWindowOn();;

			// 	// The filter in charge of effectively extracting the points is
			// 	// fed with the data.
			// 	// The filter is provided with a pointer to the vtkPolyData of interest
			// 	select_visible_points -> SetInputDataObject(this -> mainwindow -> get_asteroid() -> get_polydata());

			// 	// The dimension of the selection area is set
			// 	// The seemingly akward operation within brackets is simply a
			// 	// conversion from a std::vector<int> containing 4 items to int[4]
			// 	select_visible_points -> SetSelection(&picker -> get_dimensions_int()[0]);
			// 	select_visible_points -> Update();

			// 	// The selected points are finally obtained in the form of a vtkPolyData
			// 	this -> selected_points_polydata = select_visible_points -> GetOutput();

			// 	// If at least one vertex was selected, the selection widget will open
			// 	if (this -> selected_points_polydata -> GetNumberOfPoints() > 0) {

			// 		// A new ModifyAreaWidget is provided with the underlying shape model
			// 		// by passing a pointer to this, which already owns it

			// 		ModifyAreaWidget * modify_area_widget = new ModifyAreaWidget(this -> mainwindow,
			// 		        this);

			// 		this -> mainwindow ->  lateral_dockwidget -> setWidget(modify_area_widget);
			// 		this -> mainwindow ->  lateral_dockwidget -> show();

			// 		// An actor is added to represent the selected cells
			// 		// Another one is also created to represent the rest of the shape model
			// 		// that is not selected
			// 		modify_area_widget -> highlight_selected_cells();

			// 		// the normals of the selected cells are computed
			// 		modify_area_widget -> compute_selected_cells_average_normals();

			// 		// the id of the blob "center" is found
			// 		modify_area_widget -> find_blob_center();

			// 		// the actors owned by the mainwindow are hidden. They will be
			// 		// shown again once pc_editing widget is closed.
			// 		// hiding them after highlighting the cells
			// 		// prevents flickering
			// 		this -> mainwindow -> set_actors_visibility(false);


			// 		// Back to orient mode
			// 		this -> CurrentMode = INTERACTOR_IS_ORIENT;

			// 	}

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