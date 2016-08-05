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

	this -> table = new QTableWidget(this);

	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, this);
	transform_direction_title = new QLabel("Transform direction:", this);
	interpolation_type_title = new QLabel("Interpolation type:", this);
	transform_selection_title = new QLabel("Transform selection:", this);
	slider_title = new QLabel("Transform magnitude:", this);
	button_show_vertex_table = new QPushButton("Show vertex table", this);

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
	main_layout -> addWidget(button_show_vertex_table, Qt::AlignCenter);
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
	connect(button_show_vertex_table, SIGNAL(clicked()), this, SLOT(show_vertex_table()));

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

}

void SelectedPointWidget::highlight_selected_cells() {

	this -> selected_cells_polydata = this -> get_selected_blob_polydata();

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
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();

}



vtkSmartPointer<vtkPolyData> SelectedPointWidget::get_selected_blob_polydata() {

	// The selected facets are highlighted by means of a polydata representing them
	vtkSmartPointer<vtkPolyData> selected_cells_polydata = vtkPolyData::New();

	// The points stored by this polydata are all the vertices of the full shape model
	selected_cells_polydata -> SetPoints(this -> points_polydata -> GetPoints());


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

	//set values
	for (std::set<int>::iterator iter = cells_to_include_indices.begin();
	        iter != cells_to_include_indices.end(); ++iter ) {

		// indices of the poly's vertices to retain
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

	// The structured polydata is finally constructed
	selected_cells_polydata -> SetPolys(selected_polys_cell_array);

	return selected_cells_polydata;

}



void SelectedPointWidget::show_vertex_table() {
	Qt::WindowFlags flags = windowFlags();

	this -> table-> setWindowFlags(flags | Qt::WindowStaysOnTopHint);

	this -> table -> setShowGrid(true);
	this -> table -> setEditTriggers(QAbstractItemView::NoEditTriggers); //  prevents the user from editing the items

	//  The dimensions of the table are set
	this -> table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	this -> table -> setColumnCount(4);

	// The header row is filled up
	this -> labels << "ID" << "x" << "y" << "z";

	// The rest of the table is filled up
	this -> table -> setHorizontalHeaderLabels(labels);
	this -> populate_vertex_table();
	this -> table -> show();
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
	
	vtkSmartPointer<vtkPoints> points = vtkPoints::New();
	points ->DeepCopy(this -> points_polydata -> GetPoints());

	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	vtkSmartPointer<vtkIdTypeArray> visible_points_ids = vtkIdTypeArray::SafeDownCast(ids);

	for (int i = 0; i < visible_points_ids -> GetNumberOfTuples () ; ++i ) {
		double p[3];
		points -> GetPoint (* (visible_points_ids -> GetTuple (i)), p);
		double new_p[3];
		new_p[0] = (1 + float(pos) / 100) * p[0] ;
		new_p[1] = (1 + float(pos) / 100) * p[1] ;
		new_p[2] = (1 + float(pos) / 100) * p[2] ;
		points -> SetPoint(* (visible_points_ids -> GetTuple (i)), new_p);
	}

	this -> selected_cells_polydata -> SetPoints(points);
	this -> selected_cells_polydata -> Modified();
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
}

void SelectedPointWidget::accept() {
	*this -> mainwindow -> selection_widget_is_open = false;

	this -> points_polydata -> SetPoints(this -> selected_cells_polydata -> GetPoints());
	this -> points_polydata -> Modified();
	this -> remove_selected_points_actor();
	this -> mainwindow -> qvtkWidget -> GetRenderWindow() -> Render();
	this -> table -> clear();
	QDialog::accept();
}

void SelectedPointWidget::reject() {
	*this -> mainwindow -> selection_widget_is_open = false;
	this -> remove_selected_points_actor();
	this -> table -> clear();
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

