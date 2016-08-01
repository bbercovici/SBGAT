#include "interactor.h"
#include "InheritedPicker.h"


InteractorStyle::InteractorStyle() {
	this->SelectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
	this->SelectedActor = vtkSmartPointer<vtkActor>::New();
	this->SelectedActor->SetMapper(this->SelectedMapper);
}

void InteractorStyle::OnLeftButtonUp() {
	// Forward events
	vtkInteractorStyleRubberBandPick::OnLeftButtonUp();

	if (this -> mainwindow -> selectorActive == true) {

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
		select_visible_points -> SetRenderer(this -> mainwindow -> getRenderer());
		select_visible_points -> Update();

		// The selected points are finally obtained in the form of a vtkPolyData
		vtkSmartPointer<vtkPolyData> selected_points_polydata = select_visible_points -> GetOutput();
		vtkSmartPointer<vtkDataArray> ids = selected_points_polydata -> GetPointData() -> GetArray("ids");
		vtkSmartPointer<vtkIdTypeArray> visible_points_ids = vtkIdTypeArray::SafeDownCast(ids);

		// If at least one vertex was selected
		if (selected_points_polydata -> GetNumberOfPoints() > 0) {





			// The widget allowing the user to visualize the selected point and
			// define the transform to be applied is populated
			SelectedPointWidget * pc_editing_widget = new SelectedPointWidget(points_polydata, selected_points_polydata);

			// The widget is run as a modal window
			pc_editing_widget -> exec();

			vtkSmartPointer<vtkPoints> points = this -> points_polydata -> GetPoints();


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

		this -> mainwindow -> selectorActive = false;
	}
}

vtkStandardNewMacro(InteractorStyle);


void InteractorStyle::SetPoints(vtkSmartPointer<vtkPolyData> points_polydata) {
	this->points_polydata = points_polydata;
}

void InteractorStyle::SetMainWindow(MainWindow * mainwindow) {
	this->mainwindow = mainwindow;
}
