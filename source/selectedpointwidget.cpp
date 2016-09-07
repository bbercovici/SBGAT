#include "selectedpointwidget.h"
#include <vtkSphereSource.h>


SelectedPointWidget::SelectedPointWidget(QWidget * parent) : QDialog(parent) {

	// The different GUI elements are created
	slider = new QSlider(Qt::Horizontal, this);
	slider_value = new QLineEdit(this);
	slider_layout = new QHBoxLayout();
	main_layout = new QVBoxLayout();
	slider_holder_widget = new QWidget(this);
	transform_direction_list = new QComboBox(this);
	interpolation_type_list = new QComboBox(this);
	transform_selection_list = new QComboBox(this);
	table = new QTableWidget(this);
	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, this);
	transform_direction_title = new QLabel("Transform direction:", this);
	interpolation_type_title = new QLabel("Interpolation type:", this);
	transform_selection_title = new QLabel("Transform selection:", this);
	slider_title = new QLabel("Transform magnitude:", this);


	selected_cells_polydata = vtkPolyData::New();
	unselected_cells_polydata = vtkPolyData::New();

	selected_points = vtkPoints::New();
	new_selected_points_coordinates = vtkPoints::New();


	selected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	selected_polys_ids -> SetNumberOfComponents(1);


	unselected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	unselected_polys_ids -> SetNumberOfComponents(1);

	cell_ids = vtkIdList::New() ;

	averaged_normal_array = vtkDoubleArray::New();
	averaged_normal_array -> SetNumberOfComponents(3);
	averaged_normal_array -> SetNumberOfTuples(1);

	// The slider position and range are set
	slider -> setMinimum(-100);
	slider -> setMaximum(100);
	slider -> setValue(0);

	// The QLine edit allowing one to read/set the slider position is set
	slider_value -> setFixedWidth(30);
	slider_value -> setText(QString::number(0));
	slider_holder_widget -> setLayout(slider_layout);
	slider_layout -> addWidget(slider);
	slider_layout -> addWidget(slider_value);

	// Forces the QLine Edit widget to only accept integer values between -100 and 100
	slider_value -> setValidator( new QIntValidator(-100, 100, this) );

	// The slider is connected to the QLineEdit
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(show_new_slider_pos(int)) );
	connect(slider, SIGNAL(valueChanged(int)), this, SLOT(update_view(int)) );

	// Conversly, the QLineEdit is connected to the slider
	connect(slider_value, SIGNAL(textEdited(QString)), this, SLOT(set_new_slider_pos()) );

	// The table showing the selected vertex info is set up
	this -> table -> setShowGrid(true);
	this -> table -> setEditTriggers(QAbstractItemView::NoEditTriggers); //  prevents the user from editing the items
	this -> table -> setColumnCount(4);
	this -> labels << "ID" << "x" << "y" << "z";
	this -> table -> setHorizontalHeaderLabels(labels);

	// The main layout is set and filled up
	main_layout -> setSpacing(0);
	main_layout -> setMargin(0);
	main_layout -> addWidget(transform_direction_title);
	main_layout -> addWidget(transform_direction_list);
	main_layout -> addWidget(interpolation_type_title);
	main_layout -> addWidget(interpolation_type_list);
	main_layout -> addWidget(transform_selection_title);
	main_layout -> addWidget(transform_selection_list);
	main_layout -> addWidget(slider_title);
	main_layout -> addWidget(slider_holder_widget);
	main_layout -> addWidget(table);
	main_layout -> addWidget(button_box, Qt::AlignCenter);

	// The two drop-down lists are filled
	transform_direction_list -> insertItem(0, "Radial");
	transform_direction_list -> insertItem(1, "Average normal");
	transform_direction_list -> insertItem(2, "Point normals");
	interpolation_type_list -> insertItem(0, "Uniform (0th order)");
	interpolation_type_list -> insertItem(1, "Linear (1st order)");
	interpolation_type_list -> insertItem(2, "Parabolic (2nd order)");
	transform_selection_list -> insertItem(0, "Selected blob");
	transform_selection_list -> insertItem(1, "N closest neighbors from center");

	// The state of those drop-down lists are set
	transform_direction = TransformDirection::RADIAL;
	interpolation_type = InterpolationType::UNIFORM;
	transform_selection = TransformSelection::SELECTED;

	// Each drop down lists generates a signal notyfing the program that it was changed
	connect( transform_direction_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_transform_direction);

	connect( interpolation_type_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_interpolation_type);

	connect( transform_selection_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_transform_selection);



	// The different buttons are connected to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));

	this -> setLayout(main_layout);
	this -> mainwindow =  qobject_cast<MainWindow*>(this -> parent());
}

void SelectedPointWidget::set_data(vtkSmartPointer<InteractorStyle> interactor_style) {

	// The selected points and the full point facet/vertex shape model are made accessible to the widget

	this -> selected_points_polydata = interactor_style -> get_selected_points_polydata();
	this -> all_points_polydata = interactor_style -> get_all_points_polydata();

	// Get the polys connectivity of the full shape model. Those are not changing,
	// and can hence be set when the new shape data is loaded

	this -> polys_ids  = this -> all_points_polydata -> GetPolys () -> GetData ();

	// Likewise, the ids of the selected points are retrieved
	this -> visible_points_global_ids_from_local_index = vtkIdTypeArray::SafeDownCast(
	            this -> selected_points_polydata -> GetPointData() -> GetArray("ids"));

	// This prevents another instance of the widget to be opened
	*this -> mainwindow -> selection_widget_is_open = true;

	// The table showing the vertex info is populated and shown
	this -> table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	this -> populate_vertex_table();
	this -> table -> show();
	this -> slider -> setValue(0);


}

void SelectedPointWidget::highlight_selected_cells() {

	// the blobs representing the selected and unselected cells are created
	this -> compute_cell_blobs();

	// The selected facets/vertices are highlighted
	vtkSmartPointer<vtkPolyDataMapper> selected_cell_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	selected_cell_mapper -> SetInputData(this -> selected_cells_polydata );
	vtkSmartPointer<vtkActor> selected_cells_actor = vtkSmartPointer<vtkActor>::New();
	selected_cells_actor -> SetMapper(selected_cell_mapper);
	selected_cells_actor -> GetMapper() -> ScalarVisibilityOff();
	selected_cells_actor -> GetProperty() -> SetColor(1, 0 , 0);
	this -> mainwindow -> get_renderer() -> AddActor(selected_cells_actor);

	// The pointer to the newly created actor is saved for future use
	this -> actor_vector.push_back(selected_cells_actor);

	// The unselected facets/vertices are also added
	vtkSmartPointer<vtkPolyDataMapper> unselected_cell_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	unselected_cell_mapper -> SetInputData(this -> unselected_cells_polydata );
	vtkSmartPointer<vtkActor> unselected_cells_actor = vtkSmartPointer<vtkActor>::New();
	unselected_cells_actor -> SetMapper(unselected_cell_mapper);
	unselected_cells_actor -> GetMapper() -> ScalarVisibilityOff();
	unselected_cells_actor -> GetProperty() -> SetColor(1, 1 , 1);
	this -> mainwindow -> get_renderer() -> AddActor(unselected_cells_actor);
	this -> actor_vector.push_back(unselected_cells_actor);
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();

}


void SelectedPointWidget::compute_cell_blobs() {

	/*************************************************************************/
	// The selected facets are highlighted by means of a polydata representing them
	// The same process is done with the unselected cells
	// The points stored by those polydatas are all the vertices of the full shape model
	selected_cells_polydata -> SetPoints(this -> all_points_polydata -> GetPoints());
	unselected_cells_polydata -> SetPoints(this -> all_points_polydata -> GetPoints());

	std::set<int> cells_to_include_indices;

	for (int selected_point_index = 0;
	        selected_point_index < visible_points_global_ids_from_local_index -> GetNumberOfTuples () ;
	        ++selected_point_index ) {
		// For each visible vertex, the ids of the cells it belongs to are stored
		this -> all_points_polydata -> GetPointCells	(* visible_points_global_ids_from_local_index -> GetTuple(selected_point_index), cell_ids);

		// Those IDs are eventually transfered in a set for uniqueness
		for (int cell_id_index = 0; cell_id_index < cell_ids -> GetNumberOfIds(); ++cell_id_index) {
			cells_to_include_indices.insert(cell_ids -> GetId (cell_id_index));
		}
	}
	/*************************************************************************/

	// Those two containers must be re-created everytime compute_cell_blobs is called. Calling Initialize on them
	// somehow prevents them to be reused. Not calling Initialize will make compute_normals crash, probably
	// because of an indexing discrepancy
	vtkSmartPointer<vtkCellArray> selected_polys_cell_array = vtkCellArray::New();
	vtkSmartPointer<vtkCellArray> unselected_polys_cell_array = vtkCellArray::New();

	/*************************************************************************/
	//The following constructs the set of cells that are not included in the selection
	// This is done by simply taking the difference between two sets:
	// - all_cells_indices: set containing the indices of all cells (all integers from zero to
	// points_polydata -> GetNumberOfCells() - 1 )
	// - cells_to_include_indices: set containing the indices of the selected cells
	std::set<int> cells_not_included_indices;
	std::set<int> all_cells_indices;


	for (int i = 0; i < this -> all_points_polydata -> GetNumberOfCells(); ++i) {
		all_cells_indices.insert(i);

	}

	// cells_not_included_indices contains the indices of the cells that were not included in the selection
	std::set_difference(all_cells_indices.begin(), all_cells_indices.end(),
	                    cells_to_include_indices.begin(), cells_to_include_indices.end(),
	                    std::inserter(cells_not_included_indices, cells_not_included_indices.end()));
	/*************************************************************************/


	for (std::set<int>::iterator iter = cells_to_include_indices.begin();
	        iter != cells_to_include_indices.end(); ++iter ) {
		// indices of the poly's vertices that were selected
		int facet_vertices_ids[3];
		facet_vertices_ids[0] = * (polys_ids -> GetTuple (4 * (*iter) + 1));
		facet_vertices_ids[1] = * (polys_ids -> GetTuple (4 * (*iter) + 2));
		facet_vertices_ids[2] = * (polys_ids -> GetTuple (4 * (*iter) + 3));
		selected_polys_ids -> InsertNextValue(3);
		selected_polys_ids -> InsertNextValue(facet_vertices_ids[0]);
		selected_polys_ids -> InsertNextValue(facet_vertices_ids[1]);
		selected_polys_ids -> InsertNextValue(facet_vertices_ids[2]);
	}

	selected_polys_cell_array -> SetCells(cells_to_include_indices.size(), selected_polys_ids);


	for (std::set<int>::iterator iter = cells_not_included_indices.begin();
	        iter != cells_not_included_indices.end(); ++iter ) {
		// indices of the poly's vertices that were not selected
		int facet_vertices_ids[3];
		facet_vertices_ids[0] = * (polys_ids -> GetTuple (4 * (*iter) + 1));
		facet_vertices_ids[1] = * (polys_ids -> GetTuple (4 * (*iter) + 2));
		facet_vertices_ids[2] = * (polys_ids -> GetTuple (4 * (*iter) + 3));
		unselected_polys_ids -> InsertNextValue(3);
		unselected_polys_ids -> InsertNextValue(facet_vertices_ids[0]);
		unselected_polys_ids -> InsertNextValue(facet_vertices_ids[1]);
		unselected_polys_ids -> InsertNextValue(facet_vertices_ids[2]);
	}

	unselected_polys_cell_array -> SetCells(cells_not_included_indices.size(), unselected_polys_ids);

	// The structured polydata of the selected cells is finally constructed
	this -> selected_cells_polydata -> SetPolys(selected_polys_cell_array);

	// same for the unselected cells
	this -> unselected_cells_polydata -> SetPolys(unselected_polys_cell_array);

}



void SelectedPointWidget::populate_vertex_table() {

	// Ids of selected points
	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	for (int row = 0; row < this -> selected_points_polydata -> GetNumberOfPoints(); ++row) {

		// The ids of each selected point is added in the first column
		QTableWidgetItem * id_item = new QTableWidgetItem(tr("%1").arg(* (ids -> GetTuple (row))));

		table -> setItem(row, 0, id_item);

		// The coordinates of each selected point are added to the other columns
		double p[3];
		this -> selected_points_polydata -> GetPoint(row, p);

		QTableWidgetItem * x_item = new QTableWidgetItem(tr("%1").arg(p[0]));
		QTableWidgetItem * y_item = new QTableWidgetItem(tr("%1").arg(p[1]));
		QTableWidgetItem * z_item = new QTableWidgetItem(tr("%1").arg(p[2]));

		table -> setItem(row, 1, x_item);
		table -> setItem(row, 2, y_item);
		table -> setItem(row, 3, z_item);
	}
}

void SelectedPointWidget::find_blob_center() {
	// The location of the blob center is found by looping over each selected point
	double blob_average_position[3];
	double blob_center_position[3];

	double point_position[3];

	for (int i = 0; i < this -> selected_cells_polydata -> GetPoints() -> GetNumberOfPoints(); ++i ) {
		this -> selected_cells_polydata -> GetPoint (i, point_position);
		blob_average_position[0] = blob_average_position[0] + point_position[0];
		blob_average_position[1] = blob_average_position[1] + point_position[1];
		blob_average_position[2] = blob_average_position[2] + point_position[2];
	}

	blob_average_position[0] = blob_average_position[0] / this -> selected_cells_polydata -> GetPoints() -> GetNumberOfPoints();
	blob_average_position[1] = blob_average_position[1] / this -> selected_cells_polydata -> GetPoints() -> GetNumberOfPoints();
	blob_average_position[2] = blob_average_position[2] / this -> selected_cells_polydata -> GetPoints() -> GetNumberOfPoints();

	// blob_average_position now contains the average location of the selected points. Note that this does not correspond to
	// a point physically present in the blob. The next step thus consists in finding the selected point that is closest to this average.
	// this point physically present in the blob is thus called the "blob center"

	this -> selected_cells_polydata -> GetPoint(this -> selected_cells_polydata -> FindPoint(blob_average_position),
	        blob_center_position);

	// The coordinates are stored in two dedicated std::vector<double>
	this -> blob_center_position.clear();
	this -> blob_average_position.clear();

	this -> blob_center_position.push_back(blob_center_position[0]);
	this -> blob_center_position.push_back(blob_center_position[1])
	this -> blob_center_position.push_back(blob_center_position[2])


	this -> blob_average_position.push_back(blob_average_position[0]);
	this -> blob_average_position.push_back(blob_average_position[1]);
	this -> blob_average_position.push_back(blob_average_position[2]);


}

void SelectedPointWidget::update_view(int pos) {

	this -> selected_points -> DeepCopy(this -> all_points_polydata -> GetPoints());

	// A scaling factor is computed to reflect the retained interpolation type.
	double scaling_power;

	switch (this -> interpolation_type ) {
	case InterpolationType::UNIFORM:
		scaling_power = 0;
		break;
	case InterpolationType::LINEAR:
		scaling_power = 1;
		break;
	case InterpolationType::PARABOLIC:
		scaling_power = 2;
		break;
	}

	for (int i = 0; i < this -> visible_points_global_ids_from_local_index -> GetNumberOfTuples () ; ++i ) {
		double p[3];
		double new_p[3];
		double transform_direction[3];

		// The i-th visible point is queried
		this -> selected_points -> GetPoint (* (this -> visible_points_global_ids_from_local_index -> GetTuple (i)), p);

		switch (this -> transform_direction) {
		case TransformDirection::RADIAL:
			new_p[0] = (1 + float(pos) / 100) * p[0] ;
			new_p[1] = (1 + float(pos) / 100) * p[1] ;
			new_p[2] = (1 + float(pos) / 100) * p[2] ;
			break;

		case TransformDirection::NORMAL_POINT:

			this -> selected_cells_normals -> GetTuple(i,
			        transform_direction);

			new_p[0] =  p[0] + transform_direction[0] * float(pos) / 100;
			new_p[1] =  p[1] + transform_direction[1] * float(pos) / 100;
			new_p[2] =  p[2] + transform_direction[2] * float(pos) / 100;
			break;

		case TransformDirection::NORMAL_AVERAGED:

			transform_direction[0] = this -> averaged_normal_array -> GetValue(0);
			transform_direction[1] = this -> averaged_normal_array -> GetValue(1);
			transform_direction[2] = this -> averaged_normal_array -> GetValue(2);

			new_p[0] =  p[0] + transform_direction[0] * float(pos) / 100;
			new_p[1] =  p[1] + transform_direction[1] * float(pos) / 100;
			new_p[2] =  p[2] + transform_direction[2] * float(pos) / 100;
			break;

		}

		this -> selected_points -> SetPoint(* (this -> visible_points_global_ids_from_local_index -> GetTuple (i)), new_p);
	}


	this -> selected_cells_polydata -> SetPoints(this -> selected_points  );
	this -> selected_cells_polydata -> Modified();
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
}

void SelectedPointWidget::accept() {
	new_selected_points_coordinates -> DeepCopy(this -> selected_cells_polydata -> GetPoints());
	this -> all_points_polydata -> SetPoints(new_selected_points_coordinates);
	this -> all_points_polydata -> Modified();

	this -> reset();

	QDialog::accept();
}

void SelectedPointWidget::reject() {

	this -> reset();

	QDialog::reject();

}

void SelectedPointWidget::remove_selected_points_actor() {

	for (std::vector<vtkSmartPointer<vtkActor> >::iterator iter = this -> actor_vector.begin();
	        iter != this->actor_vector.end(); ++iter) {
		this -> mainwindow -> get_renderer() -> RemoveActor(*iter);
	}

	actor_vector.clear();
}

void SelectedPointWidget::show_new_slider_pos(int pos) {
	this -> slider_value -> setText( QString::number(pos));
}

void SelectedPointWidget::set_new_slider_pos() {
	this -> slider -> setValue(this -> slider_value -> text().toInt());

}

void SelectedPointWidget::reset() {

	*this -> mainwindow -> selection_widget_is_open = false;
	this -> remove_selected_points_actor();

	this -> table -> clear();


	// Do not call initialize on this -> selected_cells_polydata !
	this -> selected_cells_polydata -> Initialize();
	this -> unselected_cells_polydata -> Initialize();
	this -> selected_polys_ids -> Initialize();
	this -> unselected_polys_ids -> Initialize();
	this -> cell_ids -> Initialize();
	this -> selected_points -> Initialize();

	this -> mainwindow -> set_actors_visibility(true);
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();

	this -> mainwindow -> leak_tracker -> PrintCurrentLeaks();

}

float SelectedPointWidget::get_actual_memory_size() {

	float memory_used = float(
	                        this -> selected_cells_polydata -> GetActualMemorySize() +
	                        this -> unselected_cells_polydata -> GetActualMemorySize() +
	                        this -> selected_polys_ids -> GetActualMemorySize() +
	                        this -> unselected_polys_ids  -> GetActualMemorySize() +
	                        this -> selected_points -> GetActualMemorySize()
	                    ) / 1024 ;

	return memory_used;
}

void SelectedPointWidget::compute_selected_cells_normals() {

	vtkSmartPointer<vtkPolyDataNormals> selected_polydata_normals_filter = vtkPolyDataNormals::New();
	selected_polydata_normals_filter -> ComputePointNormalsOff();
	selected_polydata_normals_filter -> ComputeCellNormalsOn();
	selected_polydata_normals_filter -> SetInputData(this -> selected_cells_polydata);
	selected_polydata_normals_filter -> Update ();


	// Cell normals. There should be as many cell normals as cells in selected_cells_polydata
	this -> selected_cells_normals = selected_polydata_normals_filter -> GetOutput()
	                                 -> GetCellData() -> GetNormals();


	// *******************************************************
	// The direction of the normals is averaged
	double average_normal_direction[3];
	for (int i = 0; i < this -> selected_cells_normals -> GetNumberOfTuples(); ++i) {
		double average_normal_direction_buffer[3];

		this -> selected_cells_normals -> GetTuple(i,
		        average_normal_direction_buffer);

		average_normal_direction[0] = average_normal_direction[0] + average_normal_direction_buffer[0];
		average_normal_direction[1] = average_normal_direction[1] + average_normal_direction_buffer[1];
		average_normal_direction[2] = average_normal_direction[2] + average_normal_direction_buffer[2];

	}

	average_normal_direction[0] = average_normal_direction[0] /  this -> selected_cells_normals -> GetNumberOfTuples();
	average_normal_direction[1] = average_normal_direction[1] /  this -> selected_cells_normals -> GetNumberOfTuples();
	average_normal_direction[2] = average_normal_direction[2] /  this -> selected_cells_normals -> GetNumberOfTuples();
	// *******************************************************

	// The averaged normal is stored for later reuse
	this -> averaged_normal_array -> SetTuple(0, average_normal_direction);



}

void SelectedPointWidget::set_transform_direction(const int item_index) {
	switch (item_index) {
	case 0:
		transform_direction = TransformDirection::RADIAL;
		break;
	case 1:
		transform_direction = TransformDirection::NORMAL_AVERAGED;
		break;
	case 2:

		transform_direction = TransformDirection::NORMAL_POINT;

		break;
	default:
		std::cout << " Case not implemented in set_transform_direction. Got item_index== " << item_index << std::endl;
		break;

	}
}


void SelectedPointWidget::set_interpolation_type(const int item_index) {
	switch (item_index) {
	case 0:
		interpolation_type = InterpolationType::UNIFORM;
		break;
	case 1:

		interpolation_type = InterpolationType::LINEAR;
		break;
	case 2:

		interpolation_type = InterpolationType::PARABOLIC;
		break;
	default:
		std::cout << " Case not implemented in set_interpolation_type. Got item_index== " << item_index << std::endl;
		break;

	}
}


void SelectedPointWidget::set_transform_selection(const int item_index) {
	switch (item_index) {
	case 0:
		transform_selection = TransformSelection::SELECTED;
		break;
	case 1:

		transform_selection = TransformSelection::NCLOSEST;
		break;
	default:
		std::cout << " Case not implemented in set_transform_selection. Got item_index== " << item_index << std::endl;
		break;

	}
}