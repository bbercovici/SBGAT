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

			vtkSmartPointer<vtkSelectVisiblePoints> select_visible_points =
			    vtkSmartPointer<vtkSelectVisiblePoints>::New();

			vtkSmartPointer<vtkHardwareSelector> selector = vtkSmartPointer<vtkHardwareSelector>::New();
			selector -> SetRenderer(this -> mainwindow -> get_renderer());


			selector -> SetArea(&picker -> get_dimensions_uint()[0]);
			selector -> SetFieldAssociation(vtkDataObject::FIELD_ASSOCIATION_CELLS);
			vtkSelection* selection = selector->Select();
			std::cout << "Selection has " << selection->GetNumberOfNodes() << " nodes." << std::endl;


			// The filter renderer is set to the current renderer
			select_visible_points -> SetRenderer(this -> mainwindow -> get_renderer());
			select_visible_points -> SelectionWindowOn();;

			// The filter in charge of effectively extracting the points is
			// fed with the data.
			// The filter is provided with a pointer to the vtkPolyData of interest
			select_visible_points -> SetInputDataObject(this -> mainwindow -> get_asteroid() -> get_polydata());

			// The dimension of the selection area is set
			// The seemingly akward operation within brackets is simply a
			// conversion from a std::vector<int> containing 4 items to int[4]
			select_visible_points -> SetSelection(&picker -> get_dimensions_int()[0]);
			select_visible_points -> Update();

			// The selected points are finally obtained in the form of a vtkPolyData
			this -> selected_points_polydata = select_visible_points -> GetOutput();

			// If at least one vertex was selected, the selection widget will open
			if (this -> selected_points_polydata -> GetNumberOfPoints() > 0) {

				// A new ModifyAreaWidget is provided with the underlying shape model
				// by passing a pointer to this, which already owns it

				ModifyAreaWidget * modify_area_widget = new ModifyAreaWidget(this -> mainwindow,
				        this);

				this -> mainwindow ->  lateral_dockwidget -> setWidget(modify_area_widget);
				this -> mainwindow ->  lateral_dockwidget -> show();

				// An actor is added to represent the selected cells
				// Another one is also created to represent the rest of the shape model
				// that is not selected
				modify_area_widget -> highlight_selected_cells();

				// the normals of the selected cells are computed
				modify_area_widget -> compute_selected_cells_normals();

				// the id of the blob "center" is found
				modify_area_widget -> find_blob_center();

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