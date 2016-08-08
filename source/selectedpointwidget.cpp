#include "selectedpointwidget.h"

SelectedPointWidget::SelectedPointWidget() {

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
	transform_direction_list -> insertItem(1, "Along blob normal");
	interpolation_type_list -> insertItem(0, "Uniform (0th order)");
	interpolation_type_list -> insertItem(1, "Linear (1st order)");
	interpolation_type_list -> insertItem(2, "Parabolic (2nd order)");
	transform_selection_list -> insertItem(0, "Selected blob");
	transform_selection_list -> insertItem(1, "N closest neighbors from center");

	// The different buttons are connected to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));

	// The layout of the widget is set
	this -> setLayout(main_layout);



}


void SelectedPointWidget::set_data(vtkSmartPointer<InteractorStyle> interactor_style) {
	// The selected points and the full point facet/vertex shape model are made accessible to the widget
	this -> selected_points_polydata = interactor_style -> get_selected_points_polydata();
	this -> points_polydata = interactor_style -> get_points_polydata();
	this -> mainwindow = interactor_style -> get_mainwindow();

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

	// The selected facets are highlighted by means of a polydata representing them
	// The same process is done with the unselected cells
	// The points stored by those polydatas are all the vertices of the full shape model
	selected_cells_polydata -> SetPoints(this -> points_polydata -> GetPoints());
	unselected_cells_polydata -> SetPoints(this -> points_polydata -> GetPoints());

	// Get the polys connectivity of the full shape model
	vtkSmartPointer<vtkCellArray> polys = this -> points_polydata -> GetPolys ();
	vtkSmartPointer<vtkIdTypeArray> polys_ids = polys -> GetData ();

	// The topology of the shape is constructed (maybe should move this to load_obj?)
	this -> points_polydata -> BuildLinks();

	// Get the ids of the selected points
	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	vtkSmartPointer<vtkIdTypeArray> visible_v_ids = vtkIdTypeArray::SafeDownCast(ids);

	vtkSmartPointer<vtkIdList> cellIds = vtkIdList::New() ;
	std::set<int> cells_to_include_indices;

	for (int selected_v_index = 0;
	        selected_v_index < visible_v_ids -> GetNumberOfTuples () ;
	        ++selected_v_index ) {
		// For each visible vertex, the ids of the cells it belongs to are stored
		this -> points_polydata -> GetPointCells	(* visible_v_ids -> GetTuple(selected_v_index), cellIds);

		// Those IDs are eventually transfered in a set for uniqueness
		for (int cell_id_index = 0; cell_id_index < cellIds -> GetNumberOfIds(); ++cell_id_index) {
			cells_to_include_indices.insert(cellIds -> GetId (cell_id_index));
		}
	}

	// The selected polys are then provided to selected_cells_polydata
	vtkSmartPointer<vtkCellArray> selected_polys_cell_array = vtkCellArray::New();
	vtkSmartPointer<vtkIdTypeArray> selected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	selected_polys_ids -> SetNumberOfComponents(1);


	// The same thing is done with the unselected cells
	vtkSmartPointer<vtkCellArray> unselected_polys_cell_array = vtkCellArray::New();
	vtkSmartPointer<vtkIdTypeArray> unselected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	unselected_polys_ids -> SetNumberOfComponents(1);


	//The following constructs the set of cells that are not included in the selection
	// This is done by simply taking the difference between two sets:
	// - all_cells_indices: set containing the indices of all cells (all integers from zero to
	// points_polydata -> GetNumberOfCells() - 1 )
	// - cells_to_include_indices: set containing the indices of the selected cells

	std::set<int> cells_not_included_indices;
	std::set<int> all_cells_indices;

	for (int i = 0; i < this -> points_polydata -> GetNumberOfCells(); ++i) {
		all_cells_indices.insert(i);
	}

	// cells_not_included_indices contains the indices of the cells that were not included in the selection
	std::set_difference(all_cells_indices.begin(), all_cells_indices.end(),
	                    cells_to_include_indices.begin(), cells_to_include_indices.end(),
	                    std::inserter(cells_not_included_indices, cells_not_included_indices.end()));


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
	selected_cells_polydata -> SetPolys(selected_polys_cell_array);


	// same for the unselected cells
	unselected_cells_polydata -> SetPolys(unselected_polys_cell_array);

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

void SelectedPointWidget::update_view(int pos) {

	this -> selected_points -> DeepCopy(this -> points_polydata -> GetPoints());

	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	vtkSmartPointer<vtkIdTypeArray> visible_points_ids = vtkIdTypeArray::SafeDownCast(ids);

	for (int i = 0; i < visible_points_ids -> GetNumberOfTuples () ; ++i ) {
		double p[3];
		selected_points -> GetPoint (* (visible_points_ids -> GetTuple (i)), p);
		double new_p[3];
		new_p[0] = (1 + float(pos) / 100) * p[0] ;
		new_p[1] = (1 + float(pos) / 100) * p[1] ;
		new_p[2] = (1 + float(pos) / 100) * p[2] ;
		selected_points -> SetPoint(* (visible_points_ids -> GetTuple (i)), new_p);
	}


	this -> selected_cells_polydata -> SetPoints(selected_points);
	this -> selected_cells_polydata -> Modified();
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
}

void SelectedPointWidget::accept() {
	*this -> mainwindow -> selection_widget_is_open = false;
	this -> points_polydata -> SetPoints(this -> selected_cells_polydata -> GetPoints());
	this -> points_polydata -> Modified();
	this -> remove_selected_points_actor();
	this -> mainwindow -> set_actors_visibility(true);
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();

	this -> table -> clear();

	QDialog::accept();
}

void SelectedPointWidget::reject() {
	*this -> mainwindow -> selection_widget_is_open = false;
	this -> remove_selected_points_actor();
	this -> table -> clear();

	this -> mainwindow -> set_actors_visibility(true);
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
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

