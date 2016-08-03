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
	this -> widget_is_open = new bool();
	*this -> widget_is_open = false;
}


void InteractorStyle::OnLeftButtonUp() {
	// Forward events

	vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

	if (this -> CurrentMode == INTERACTOR_IS_SELECT && *this -> widget_is_open == false) {


		InheritedPicker * picker ;
		picker = static_cast< InheritedPicker * >(this -> Interactor -> GetPicker());


		// The filter in charge of effectively extracting the points is
		// created.
		vtkSmartPointer<vtkSelectVisiblePoints> select_visible_points =
		    vtkSmartPointer<vtkSelectVisiblePoints>::New();

		select_visible_points -> SelectionWindowOn();

		// The filter is provided with a pointer to the vtkPolyData of interest
		select_visible_points -> SetInputDataObject(this -> points_polydata);

		// The dimension of the selection area is set
		// The seemingly akward operation within brackets is simply a
		// conversion from a std::vector<int> containing 4 items to int[4]
		select_visible_points -> SetSelection(&picker -> get_dimensions()[0]);

		// The filter renderer is set to the current renderer
		select_visible_points -> SetRenderer(this -> mainwindow -> get_renderer());
		select_visible_points -> Update();

		


		// The selected points are finally obtained in the form of a vtkPolyData
		this -> selected_points_polydata = select_visible_points -> GetOutput();
		// this -> selected_cells_polydata = select_visible_cells -> GetOutput();

		// If at least one vertex was selected, the selection widget will open
		if (this -> selected_points_polydata -> GetNumberOfPoints() > 0) {


			// The widget allowing the user to visualize the selected point and
			// define the transform to be applied is created
			
			SelectedPointWidget * pc_editing_widget = new SelectedPointWidget(this);

			// The widget is run as a non-modal window
			pc_editing_widget -> show();

			// Back to orient mode
			this -> CurrentMode = INTERACTOR_IS_ORIENT;

		}



	}

}

vtkStandardNewMacro(InteractorStyle);

void InteractorStyle::transform_points() {
	vtkSmartPointer<vtkPoints> points = this -> points_polydata -> GetPoints();
	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	vtkSmartPointer<vtkIdTypeArray> visible_points_ids = vtkIdTypeArray::SafeDownCast(ids);

	for (int i = 0; i < visible_points_ids -> GetNumberOfTuples () ; ++i ) {
		double p[3];
		points -> GetPoint (* (visible_points_ids -> GetTuple (i)), p);
		double new_p[3];
		new_p[0] = 0.75 * p[0] ;
		new_p[1] = 0.75 * p[1] ;
		new_p[2] = 0.75 * p[2] ;
		points -> SetPoint(* (visible_points_ids -> GetTuple (i)), new_p);

	}
	this -> points_polydata -> SetPoints(points);
	this -> points_polydata -> Modified();
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
}

void InteractorStyle::set_points(vtkSmartPointer<vtkPolyData> points_polydata) {
	this -> points_polydata = points_polydata;
}

void InteractorStyle::set_mainwindow(MainWindow * mainwindow) {
	this->mainwindow = mainwindow;
}

void InteractorStyle::set_current_mode(const int mode) {
	this -> CurrentMode = mode;
}

InteractorStyle::~InteractorStyle() {
	delete (this -> widget_is_open);
}

MainWindow * InteractorStyle::get_mainwindow(){
	return this -> mainwindow;
}


vtkSmartPointer<vtkPolyData> InteractorStyle::get_selected_points_polydata(){
	return this -> selected_points_polydata;
}
vtkSmartPointer<vtkPolyData> InteractorStyle::get_points_polydata(){
	return this -> points_polydata;
}

bool * InteractorStyle::get_widget_is_open(){
	return this -> widget_is_open;
}
