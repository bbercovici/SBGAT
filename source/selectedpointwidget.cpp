#include "selectedpointwidget.h"

SelectedPointWidget::SelectedPointWidget(vtkSmartPointer<InteractorStyle> interactor_style) {

	// The different GUI elements are created
	this -> setAttribute(Qt::WA_DeleteOnClose);

	layout = new QHBoxLayout();
	list_holder_layout = new QVBoxLayout();
	list_holder_widget = new QWidget();
	transform_direction_list = new QComboBox();
	interpolation_type_list = new QComboBox();
	transform_selection_list = new QComboBox();

	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Vertical);
	transform_direction_title = new QLabel("Transform direction:");
	interpolation_type_title = new QLabel("Interpolation type:");
	transform_selection_title = new QLabel("Transform selection:");

	button_show_vertex_table = new QPushButton("Show vertex table");

	// The drop-down lists are added to the corresponding layout
	list_holder_widget -> setLayout(list_holder_layout);
	list_holder_layout -> setSpacing(0);
	list_holder_layout -> setMargin(0);
	list_holder_layout -> addWidget(transform_direction_title, Qt::AlignTop);
	list_holder_layout -> addWidget(transform_direction_list, Qt::AlignTop);
	list_holder_layout -> addWidget(interpolation_type_title, Qt::AlignTop);
	list_holder_layout -> addWidget(interpolation_type_list, Qt::AlignTop);
	list_holder_layout -> addWidget(transform_selection_title, Qt::AlignTop);
	list_holder_layout -> addWidget(transform_selection_list, Qt::AlignTop);



	// The two drop-down lists are filled
	transform_direction_list -> insertItem(0, "Radial");
	transform_direction_list -> insertItem(1, "Along blob normal");
	interpolation_type_list -> insertItem(0, "Uniform (0th order)");
	interpolation_type_list -> insertItem(1, "Linear (1st order)");
	interpolation_type_list -> insertItem(2, "Parabolic (2nd order)");
	transform_selection_list -> insertItem(0, "Selected blob");
	transform_selection_list -> insertItem(1, "N closest neighbors from center");

	// The selected points and the full point facet/vertex shape model are made accessible to the widget
	this -> selected_points_polydata = interactor_style -> get_selected_points_polydata();
	this -> points_polydata = interactor_style -> get_points_polydata();

	// The different widgets are added to the outer layout
	layout -> addWidget(list_holder_widget);
	layout -> addWidget(button_box, Qt::AlignCenter);
	layout -> addWidget(button_show_vertex_table, Qt::AlignCenter);
	this -> setLayout(layout);

	// The different buttons are connected to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));
	connect(button_show_vertex_table, SIGNAL(clicked()), this, SLOT(show_vertex_table()));


	// This forces the widget to stay on top of all others
	Qt::WindowFlags flags = windowFlags();
	this -> setWindowFlags(flags | Qt::WindowStaysOnTopHint);

	// This prevents another instance of the widget to be opened
	this -> widget_is_open = interactor_style -> get_widget_is_open();
	*this -> widget_is_open = true;

	// This links the selection widget to the interactor
	this -> interactor_style = interactor_style;

	// The selected facets/vertices are highlighted
	vtkSmartPointer<vtkPolyDataMapper> selected_cell_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	selected_cell_mapper -> SetInputData(this -> get_selected_blob_polydata() );
	vtkSmartPointer<vtkActor> selected_cells_actor = vtkSmartPointer<vtkActor>::New();

	selected_cells_actor -> SetMapper(selected_cell_mapper);
	selected_cells_actor -> GetMapper() -> ScalarVisibilityOff();
	selected_cells_actor -> GetProperty() -> SetColor(1, 0 , 0);

	this -> interactor_style -> get_mainwindow() -> get_renderer() -> AddActor(selected_cells_actor);

}





vtkSmartPointer<vtkPolyData> SelectedPointWidget::get_selected_blob_polydata() {

	// The selected facets are highlighted by means of a polydata representing them
	vtkSmartPointer<vtkPolyData> selected_cells_polydata = vtkPolyData::New();

	// The points stored by this polydata are all the vertices of the full shape model
	selected_cells_polydata -> SetPoints(this -> interactor_style -> get_points_polydata() -> GetPoints());


	// Get the polys connectivity of the full shape model
	vtkSmartPointer<vtkCellArray> polys = interactor_style -> get_points_polydata() -> GetPolys ();
	vtkSmartPointer<vtkIdTypeArray> polys_ids = polys -> GetData ();

	// The topology of the shape is constructed (maybe should move this to load_obj?)
	this -> interactor_style -> get_points_polydata() -> BuildLinks();


	// Get the ids of the selected points
	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	vtkSmartPointer<vtkIdTypeArray> visible_v_ids = vtkIdTypeArray::SafeDownCast(ids);

	vtkSmartPointer<vtkIdList> cellIds = vtkIdList::New() ;
	std::set<int> cells_to_include_indices;

	for (int selected_v_index = 0;
	        selected_v_index < visible_v_ids -> GetNumberOfTuples () ;
	        ++selected_v_index ) {
		// For each visible vertex, the ids of the cells it belongs to are stored
		this -> interactor_style ->
		get_points_polydata() -> GetPointCells	(* visible_v_ids -> GetTuple(selected_v_index), cellIds);

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

	this -> table = new QTableWidget(this);
	this -> table-> setWindowFlags(flags | Qt::WindowStaysOnTopHint);

	this -> table -> setShowGrid(true);
	this -> table -> setEditTriggers(QAbstractItemView::NoEditTriggers); //  prevents the user from editing the items

	//  The dimensions of the table are set
	this -> table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	this -> table -> setColumnCount(4);

	this -> labels << "ID" << "x" << "y" << "z";
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

void SelectedPointWidget::accept() {
	(*this -> widget_is_open) = false;
	this -> interactor_style -> transform_points();
	this -> remove_selected_points_actor();
	this -> interactor_style -> get_mainwindow()-> qvtkWidget -> GetRenderWindow() -> Render();

	QDialog::accept();


}

void SelectedPointWidget::reject() {
	(*this -> widget_is_open) = false;
	this -> remove_selected_points_actor();
	this -> interactor_style -> get_mainwindow()-> qvtkWidget -> GetRenderWindow() -> Render();

	QDialog::reject();

}

void SelectedPointWidget::remove_selected_points_actor() {
	vtkSmartPointer<vtkActorCollection> actor_coll = this -> interactor_style -> get_mainwindow()
	        -> get_renderer() -> GetActors();
	this -> interactor_style -> get_mainwindow()
	-> get_renderer() -> RemoveActor(actor_coll -> GetLastActor());
}
