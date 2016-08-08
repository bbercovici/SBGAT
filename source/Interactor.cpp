#include "interactor.h"
#include "InheritedPicker.h"
#include "utilities.h"

// shortcut to interactor modes
#define INTERACTOR_IS_ORIENT 0
#define INTERACTOR_IS_SELECT 1


InteractorStyle::InteractorStyle() {
	this -> SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	this -> SelectedActor = vtkSmartPointer<vtkActor>::New();
	this -> SelectedActor -> SetMapper(this -> SelectedMapper);
	this -> CurrentMode = INTERACTOR_IS_ORIENT;
	this -> mainwindow = NULL;
	select_visible_points =
	    vtkSmartPointer<vtkSelectVisiblePoints>::New();
	select_visible_points -> SelectionWindowOn();
}


void InteractorStyle::OnLeftButtonUp() {
	// Forward events
	vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

	if (this -> mainwindow != NULL) {
		// If the interactor is in selection mode
		// AND if the selection widget is not already open
		if (this -> CurrentMode == INTERACTOR_IS_SELECT
		        &&  *this -> mainwindow -> selection_widget_is_open == false) {

			InheritedPicker * picker ;
			picker = static_cast< InheritedPicker * >(this -> Interactor -> GetPicker());

			// The filter in charge of effectively extracting the points is
			// fed with the data.
			// The filter is provided with a pointer to the vtkPolyData of interest
			select_visible_points -> SetInputDataObject(this -> points_polydata);

			// // The dimension of the selection area is set
			// // The seemingly akward operation within brackets is simply a
			// // conversion from a std::vector<int> containing 4 items to int[4]
			select_visible_points -> SetSelection(&picker -> get_dimensions()[0]);
			select_visible_points -> Update();

			// // The selected points are finally obtained in the form of a vtkPolyData
			this -> selected_points_polydata = select_visible_points -> GetOutput();

			// // If at least one vertex was selected, the selection widget will open
			if (this -> selected_points_polydata -> GetNumberOfPoints() > 0) {

				// 	// the widget is provided with the underlying shape model
				// 	// by passing a pointer to this, which already owns it
				this -> mainwindow -> pc_editing_widget -> set_data(this);

				// An actor is added to represent the selected cells
				// Another one is also create to represent the rest of the shape model
				// that is not selected
				this -> mainwindow -> pc_editing_widget -> highlight_selected_cells();

				// the actors owned by the mainwindow are hidden. They will be
				// shown again once pc_editing widget is closed.
				// hiding them after highlighting the cells
				// prevents flickering
				this -> mainwindow -> set_actors_visibility(false);

				// The widget is run as a non-modal window
				this -> mainwindow -> pc_editing_widget -> show();

				// Back to orient mode
				this -> CurrentMode = INTERACTOR_IS_ORIENT;

			}

		}

	}

}

vtkStandardNewMacro(InteractorStyle);


void InteractorStyle::set_points(vtkSmartPointer<vtkPolyData> points_polydata) {
	this -> points_polydata = points_polydata;
}

void InteractorStyle::set_mainwindow(MainWindow * mainwindow) {
	this -> mainwindow = mainwindow;
	// The filter renderer is set to the current renderer
	select_visible_points -> SetRenderer(this -> mainwindow -> get_renderer());

}

void InteractorStyle::set_current_mode(const int mode) {
	this -> CurrentMode = mode;
}


MainWindow * InteractorStyle::get_mainwindow() {
	return this -> mainwindow;
}


vtkSmartPointer<vtkPolyData> InteractorStyle::get_selected_points_polydata() {
	return this -> selected_points_polydata;
}
vtkSmartPointer<vtkPolyData> InteractorStyle::get_points_polydata() {
	return this -> points_polydata;
}
