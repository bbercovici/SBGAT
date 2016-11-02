#include "SelectedPointWidget.hpp"


SelectedPointWidget::SelectedPointWidget(Mainwindow * parent,
        InteractorStyle * interactor_style) : QDialog(parent) {

	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;


	// The different GUI elements are created
	this -> slider_magnitude = new QSlider(Qt::Horizontal, this);
	this -> slider_magnitude_value = new QLineEdit(this);
	this -> slider_magnitude_layout = new QHBoxLayout();

	this -> slider_neighbors = new QSlider(Qt::Horizontal, this);
	this -> slider_neighbors_value = new QLineEdit(this);
	this -> slider_neighbors_layout = new QHBoxLayout();


	this -> main_layout = new QVBoxLayout();

	this -> slider_magnitude_holder_widget = new QWidget(this);
	this -> slider_neighbors_holder_widget = new QWidget(this);

	this -> transform_direction_list = new QComboBox(this);
	this -> interpolation_type_list = new QComboBox(this);
	this -> transform_selection_list = new QComboBox(this);
	this -> table = new QTableWidget(this);
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Horizontal, this);
	this -> transform_direction_title = new QLabel("Transform direction:", this);
	this -> interpolation_type_title = new QLabel("Interpolation type:", this);
	this -> transform_selection_title = new QLabel("Transform selection:", this);
	this -> slider_magnitude_title = new QLabel("Transform magnitude:", this);
	this -> slider_neighbors_title = new QLabel("Number of neighbors:", this);

	this -> point_locator = vtkSmartPointer<vtkKdTreePointLocator>::New();

	this -> active_selected_points_polydata = vtkPolyData::New();
	this -> active_selected_points_polydata -> SetPoints(vtkPoints::New());

	this -> selected_cells_polydata = vtkPolyData::New();
	this -> unselected_cells_polydata = vtkPolyData::New();

	this -> selected_points = vtkPoints::New();
	this -> new_selected_points_coordinates = vtkPoints::New();

	this -> selected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	this -> selected_polys_ids -> SetNumberOfComponents(1);


	this -> unselected_polys_ids =
	    vtkSmartPointer<vtkIdTypeArray>::New();
	this -> unselected_polys_ids -> SetNumberOfComponents(1);

	this -> cell_ids = vtkIdList::New() ;
	this -> N_closest_vertices_indices = vtkIdList::New() ;

	this -> averaged_normal_array = vtkDoubleArray::New();
	this -> averaged_normal_array -> SetNumberOfComponents(3);
	this -> averaged_normal_array -> SetNumberOfTuples(1);

	// The slider_magnitude position and range are set
	this -> slider_magnitude -> setMinimum(-100);
	this -> slider_magnitude -> setMaximum(100);
	this -> slider_magnitude -> setValue(0);

	// The slider_neighbors position and range are set
	this -> slider_neighbors -> setMinimum(1);
	this -> slider_neighbors -> setMaximum(10);
	this -> slider_neighbors -> setValue(10);


	// The QLine edit allowing one to read/set the slider_magnitude position is set
	this -> slider_magnitude_value -> setFixedWidth(30);
	this -> slider_magnitude_value -> setText(QString::number(0));
	this -> slider_magnitude_holder_widget -> setLayout(slider_magnitude_layout);
	this -> slider_magnitude_layout -> addWidget(slider_magnitude);
	this -> slider_magnitude_layout -> addWidget(slider_magnitude_value);

	// The QLine edit allowing one to read/set the slider_neighbors position is set
	this -> slider_neighbors_value -> setFixedWidth(30);
	this -> slider_neighbors_value -> setText(QString::number(1));
	this -> slider_neighbors_holder_widget -> setLayout(slider_neighbors_layout);
	this -> slider_neighbors_layout -> addWidget(slider_neighbors);
	this -> slider_neighbors_layout -> addWidget(slider_neighbors_value);

	// Forces the QLine Edits widget to only accept integer values between the specified bounds
	this -> slider_magnitude_value -> setValidator( new QIntValidator(-100, 100, this) );
	this -> slider_neighbors_value -> setValidator( new QIntValidator(1, 100, this) );

	// The slider_magnitude is connected to the QLineEdit
	connect(this -> slider_magnitude, SIGNAL(valueChanged(int)), this, SLOT(show_new_slider_magnitude_pos(int)) );
	connect(this -> slider_magnitude, SIGNAL(valueChanged(int)), this, SLOT(update_view(int)) );

	connect(this -> slider_neighbors, SIGNAL(valueChanged(int)), this, SLOT(show_new_slider_neighbors_pos(int)) );
	connect(this -> slider_neighbors, SIGNAL(valueChanged(int)), this, SLOT(find_N_neighbors_indices_and_update_view(int)) );

	// Conversly, the QLineEdit is connected to the slider_magnitude
	connect(this -> slider_magnitude_value, SIGNAL(textEdited(QString)), this, SLOT(set_new_slider_magnitude_pos()) );
	connect(this -> slider_neighbors_value, SIGNAL(textEdited(QString)), this, SLOT(set_new_slider_neighbors_pos()) );

	// The table showing the selected vertex info is set up
	this -> table -> setShowGrid(true);
	//  prevents the user from editing the items
	this -> table -> setEditTriggers(QAbstractItemView::NoEditTriggers);
	this -> table -> setColumnCount(4);
	this -> labels << "ID" << "x" << "y" << "z";
	this -> table -> setHorizontalHeaderLabels(labels);

	// The main layout is set and filled up
	this -> main_layout -> setSpacing(0);
	this -> main_layout -> setMargin(0);
	this -> main_layout -> addWidget(this -> transform_direction_title);
	this -> main_layout -> addWidget(this -> transform_direction_list);
	this -> main_layout -> addWidget(this -> interpolation_type_title);
	this -> main_layout -> addWidget(this -> interpolation_type_list);
	this -> main_layout -> addWidget(this -> transform_selection_title);
	this -> main_layout -> addWidget(this -> transform_selection_list);
	this -> main_layout -> addWidget(this -> slider_magnitude_title);
	this -> main_layout -> addWidget(this -> slider_magnitude_holder_widget);
	this -> main_layout -> addWidget(this -> slider_neighbors_title);
	this -> main_layout -> addWidget(this -> slider_neighbors_holder_widget);

	// The neighbors level widgets are hidden
	this -> slider_neighbors_title -> setVisible(false);
	this -> slider_neighbors_holder_widget -> setVisible(false);

	this -> main_layout -> addWidget(this -> table);
	this -> main_layout -> addWidget(button_box, Qt::AlignCenter);

	// The two drop-down lists are filled
	this -> transform_direction_list -> insertItem(0, "Radial");
	this -> transform_direction_list -> insertItem(1, "Average normal");
	this -> transform_direction_list -> insertItem(2, "Point normals");
	this -> interpolation_type_list -> insertItem(0, "Uniform (0th order)");
	this -> interpolation_type_list -> insertItem(1, "Linear (1st order)");
	this -> interpolation_type_list -> insertItem(2, "Parabolic (2nd order)");
	this -> transform_selection_list -> insertItem(0, "Selected blob");
	this -> transform_selection_list -> insertItem(1, "N closest neighbors from center");

	// The state of those drop-down lists are set
	this -> transform_direction = TransformDirection::RADIAL;
	this -> interpolation_type = InterpolationType::UNIFORM;
	this -> transform_selection = TransformSelection::SELECTED;

	// Each drop down lists generates a signal notyfing the program that it was changed
	connect( transform_direction_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_transform_direction);

	connect( interpolation_type_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_interpolation_type);

	connect( transform_selection_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::set_transform_selection);

	connect( transform_direction_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::update_view_unchanged_slider_magnitude);

	connect( interpolation_type_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &SelectedPointWidget::update_view_unchanged_slider_magnitude);

	// The different buttons are connected to the corresponding slots
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> button_box, SIGNAL(rejected()), this, SLOT(reject()));

	this -> setLayout(this -> main_layout);

	this -> set_data(interactor_style);
}

void SelectedPointWidget::set_data(InteractorStyle * interactor_style) {

	// The selected points and the full point facet/vertex shape model are made accessible to the widget
	this -> selected_points_polydata = interactor_style -> get_selected_points_polydata();
	this -> all_points_polydata = interactor_style -> get_all_points_polydata();

	// Get the polys connectivity of the full shape model. Those are not changing,
	// and can hence be set when the new shape data is loaded

	this -> polys_ids  = this -> all_points_polydata -> GetPolys () -> GetData ();

	// Likewise, the ids of the selected points are retrieved
	this -> visible_points_global_ids_from_local_index = vtkIdTypeArray::SafeDownCast(
	            this -> selected_points_polydata -> GetPointData() -> GetArray("ids"));

	// The table showing the vertex info is populated and shown
	this -> table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	this -> populate_vertex_table();
	this -> table -> show();
	this -> slider_magnitude -> setValue(0);

	// The point locator is cleaned up and provided with the selected blob
	this -> point_locator -> SetDataSet(this -> selected_points_polydata);
	this -> point_locator -> BuildLocator();

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
	this -> parent -> get_renderer() -> AddActor(selected_cells_actor);

	// The pointer to the newly created actor is saved for future use
	this -> actor_vector.push_back(selected_cells_actor);

	// The unselected facets/vertices are also added
	vtkSmartPointer<vtkPolyDataMapper> unselected_cell_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	unselected_cell_mapper -> SetInputData(this -> unselected_cells_polydata );
	vtkSmartPointer<vtkActor> unselected_cells_actor = vtkSmartPointer<vtkActor>::New();
	unselected_cells_actor -> SetMapper(unselected_cell_mapper);
	unselected_cells_actor -> GetMapper() -> ScalarVisibilityOff();
	unselected_cells_actor -> GetProperty() -> SetColor(1, 1 , 1);
	this -> parent -> get_renderer() -> AddActor(unselected_cells_actor);
	this -> actor_vector.push_back(unselected_cells_actor);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

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
	double blob_average_position[3]; // selected blob average coordinates
	double blob_center_position[3]; // coordinates of vertex that's closest to the blob average
	double point_position[3]; // local variable storing the coordinates of a vertex

	for (int i = 0; i < this -> selected_points_polydata -> GetPoints() -> GetNumberOfPoints(); ++i ) {
		this -> selected_points_polydata -> GetPoint (i, point_position);
		blob_average_position[0] = blob_average_position[0] + point_position[0];
		blob_average_position[1] = blob_average_position[1] + point_position[1];
		blob_average_position[2] = blob_average_position[2] + point_position[2];
	}

	blob_average_position[0] = blob_average_position[0] / this -> selected_points_polydata -> GetPoints() -> GetNumberOfPoints();
	blob_average_position[1] = blob_average_position[1] / this -> selected_points_polydata -> GetPoints() -> GetNumberOfPoints();
	blob_average_position[2] = blob_average_position[2] / this -> selected_points_polydata -> GetPoints() -> GetNumberOfPoints();

	// blob_average_position now contains the average location of the selected points. Note that this does not correspond to
	// a point physically present in the blob. The next step thus consists in finding the selected point that is closest to this average.
	// this point physically present in the blob is thus called the "blob center"

	this -> blob_center_id = this -> selected_points_polydata -> FindPoint(blob_average_position);
	this -> selected_points_polydata -> GetPoint(this -> blob_center_id,
	        blob_center_position);

	// The computed locations are stored in arma::vec
	this -> blob_center_position = {blob_center_position[0], blob_center_position[1], blob_center_position[2] };
	this -> blob_average_position = {blob_average_position[0], blob_average_position[1], blob_average_position[2] };

}

void SelectedPointWidget::update_view(int pos) {

	// This copy ensures that the user can modify a copy of the shape model and not the shape model itself.
	this -> selected_points -> DeepCopy(this -> all_points_polydata -> GetPoints());

	// Each selected point is looped over and transformed accordingly
	for (int i = 0; i < this -> visible_points_global_ids_from_local_index -> GetNumberOfTuples () ; ++i ) {

		// local variables used to communicate with VTK routines
		double p[3];
		double new_p[3];
		double normal_direction[3];

		// interpolating factor used to locally smooth the transform
		double interpolating_factor;

		// Arma::vecs
		arma::vec u; // normalized transform direction
		arma::vec new_p_vec; // new point location
		arma::vec old_p_vec; // old point location

		// The i-th visible point is queried
		this -> selected_points -> GetPoint (* (this -> visible_points_global_ids_from_local_index -> GetTuple (i)), p);

		// its position is stored in a arma::vec to enable further manipulations
		old_p_vec = {p[0], p[1], p[2]};

		switch (this -> transform_direction) {
		case TransformDirection::RADIAL :
			u = arma::normalise(old_p_vec);
			break;

		case TransformDirection::NORMAL_POINT :
			this -> selected_cells_normals -> GetTuple(i,
			        normal_direction);

			u = {normal_direction[0], normal_direction[1], normal_direction[2] };
			u = arma::normalise(u);
			break;

		case TransformDirection::NORMAL_AVERAGED :
			normal_direction[0] = this -> averaged_normal_array -> GetValue(0);
			normal_direction[1] = this -> averaged_normal_array -> GetValue(1);
			normal_direction[2] = this -> averaged_normal_array -> GetValue(2);

			u = {normal_direction[0], normal_direction[1], normal_direction[2] };
			u = arma::normalise(u);

			break;

		}


		// A characteristic length is extracted from the polydata's dimensions
		double normalizing_constant;
		switch (this -> transform_selection) {
		case TransformSelection::SELECTED:
			normalizing_constant =  this -> selected_points_polydata -> GetLength() * 0.3;
			break;

		case TransformSelection::NCLOSEST:
			normalizing_constant =  this -> active_selected_points_polydata -> GetLength() * 0.3;
			break;
		}


		// An interpolation factor is set depending on
		// the current choice of interpolating type. the 0.3 factor appearing in the denominator
		// appears to be a good compromise
		switch (this -> interpolation_type ) {
		case InterpolationType::UNIFORM:
			interpolating_factor = 1;
			break;

		case InterpolationType::LINEAR:
			interpolating_factor = std::max(0., 1 - arma::norm(this -> blob_center_position - old_p_vec)
			                                / normalizing_constant  );

			break;
		case InterpolationType::PARABOLIC:
			interpolating_factor = std::max(0., 1 - std::pow(arma::norm(this -> blob_center_position - old_p_vec)
			                                / normalizing_constant , 2));
			break;
		}

		// The position of the transform vertex is computed
		new_p_vec = old_p_vec + u * float(pos) / 100 *  interpolating_factor;

		new_p[0] = new_p_vec(0);
		new_p[1] = new_p_vec(1);
		new_p[2] = new_p_vec(2);

		this -> selected_points -> SetPoint(* (this -> visible_points_global_ids_from_local_index -> GetTuple (i)), new_p);

	}


	this -> selected_cells_polydata -> SetPoints(this -> selected_points  );
	this -> selected_cells_polydata -> Modified();
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
}

void SelectedPointWidget::accept() {
	this -> new_selected_points_coordinates -> DeepCopy(this -> selected_cells_polydata -> GetPoints());
	this -> all_points_polydata -> SetPoints(this -> new_selected_points_coordinates);
	this -> all_points_polydata -> Modified();
	this -> close();
}

void SelectedPointWidget::reject() {
	this -> close();
}

void SelectedPointWidget::remove_selected_points_actor() {

	for (std::vector<vtkSmartPointer<vtkActor> >::iterator iter = this -> actor_vector.begin();
	        iter != this -> actor_vector.end(); ++iter) {
		this -> parent -> get_renderer() -> RemoveActor(*iter);
	}

	this ->  actor_vector.clear();
}


void SelectedPointWidget::show_new_slider_neighbors_pos(int pos) {
	this -> slider_neighbors_value -> setText( QString::number(pos));
}

void SelectedPointWidget::set_new_slider_neighbors_pos() {
	this -> slider_neighbors -> setValue(this -> slider_neighbors_value -> text().toInt());

}

void SelectedPointWidget::show_new_slider_magnitude_pos(int pos) {
	this -> slider_magnitude_value -> setText( QString::number(pos));
}

void SelectedPointWidget::set_new_slider_magnitude_pos() {
	this -> slider_magnitude -> setValue(this -> slider_magnitude_value -> text().toInt());

}

void SelectedPointWidget::close() {

	this -> remove_selected_points_actor();

	this -> parent -> lateral_dockwidget -> hide();
	this -> parent -> set_action_status(true, this -> parent -> selectPointAct);

	this -> parent -> set_actors_visibility(true);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

	QDialog::close();


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

	// The averaged normal is stored for later reuse
	this -> averaged_normal_array -> SetTuple(0, average_normal_direction);



}

void SelectedPointWidget::find_N_neighbors_indices(const int N) {
	double center_point[3];
	center_point[0] = this -> blob_center_position(0);
	center_point[1] = this -> blob_center_position(1);
	center_point[2] = this -> blob_center_position(2);

	this -> point_locator -> FindClosestNPoints(N, center_point, this -> N_closest_vertices_indices);
	this -> active_selected_points_polydata -> GetPoints() -> Initialize();
	this -> selected_points_polydata -> GetPoints() -> GetPoints(this -> N_closest_vertices_indices, this -> active_selected_points_polydata -> GetPoints());
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
		this -> slider_neighbors_title -> setVisible(false);
		this -> slider_neighbors_holder_widget -> setVisible(false);
		this -> slider_neighbors -> setValue(this -> selected_points_polydata -> GetNumberOfPoints());
		// The line above is necessary because it will call the slot updating the view

		break;
	case 1:
		transform_selection = TransformSelection::NCLOSEST;
		this -> slider_neighbors_title -> setVisible(true);
		this -> slider_neighbors_holder_widget -> setVisible(true);
		this -> slider_neighbors -> setMaximum(this -> selected_points_polydata -> GetNumberOfPoints());
		this -> slider_neighbors -> setValue(this -> selected_points_polydata -> GetNumberOfPoints());
		break;
	default:
		std::cout << " Case not implemented in set_transform_selection. Got item_index== " << item_index << std::endl;
		break;

	}
}

void SelectedPointWidget::update_view_unchanged_slider_magnitude() {
	int pos = this -> slider_magnitude -> value();
	this -> update_view(pos);
}

void SelectedPointWidget::find_N_neighbors_indices_and_update_view(const int N) {
	this -> find_N_neighbors_indices(N);
	this -> update_view_unchanged_slider_magnitude();
}