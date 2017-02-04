#include "ModifyAreaWidget.hpp"


ModifyAreaWidget::ModifyAreaWidget(Mainwindow * parent,
                                   vtkSmartPointer<vtkPolyData> selected_polydata,
                                   vtkSmartPointer<vtkPolyData> unselected_polydata,
                                   vtkSmartPointer<vtkActor> selected_actor,
                                   vtkSmartPointer<vtkActor> unselected_actor,
                                   vtkSmartPointer<vtkIdList> boundary_vertex_ids_list) : QDialog(parent) {

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
	         this, &ModifyAreaWidget::set_transform_direction);

	connect( interpolation_type_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &ModifyAreaWidget::set_interpolation_type);

	connect( transform_selection_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &ModifyAreaWidget::set_transform_selection);

	connect( transform_direction_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &ModifyAreaWidget::update_view_unchanged_slider_magnitude);

	connect( interpolation_type_list, static_cast<void(QComboBox::*)(int)>(&QComboBox::currentIndexChanged),
	         this, &ModifyAreaWidget::update_view_unchanged_slider_magnitude);

	// The different buttons are connected to the corresponding slots
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> button_box, SIGNAL(rejected()), this, SLOT(reject()));

	this -> setLayout(this -> main_layout);

	this -> boundary_vertex_ids_list = boundary_vertex_ids_list;
	this -> selected_polydata = selected_polydata;
	this -> unselected_polydata = unselected_polydata;

	this -> selected_polydata_original = vtkSmartPointer<vtkPolyData>::New();
	this -> selected_polydata_original -> DeepCopy(this -> selected_polydata);

	this -> selected_actor = selected_actor;
	this -> unselected_actor = unselected_actor;

	// The table containing the global ids of the selected vertices is populated
	this -> set_data();

	// The average normal of the selected blob is computed
	this -> compute_selected_cells_average_normals();

	// the id of the blob "center" is found
	this -> find_blob_center();

}

void ModifyAreaWidget::set_data() {

	// Likewise, the ids of the selected points are retrieved
	this -> selected_vertices_global_ids_from_local_ids = vtkIdTypeArray::SafeDownCast(
	            this -> selected_polydata -> GetPointData() -> GetArray("ids"));

	// The table showing the vertex info is populated and shown
	this -> table -> setRowCount(this -> selected_polydata -> GetNumberOfPoints());
	this -> populate_vertex_table();
	this -> table -> show();
	this -> slider_magnitude -> setValue(0);

	// The point locator is cleaned up and provided with the selected blob
	this -> point_locator -> SetDataSet(this -> selected_polydata);
	this -> point_locator -> BuildLocator();

}

void ModifyAreaWidget::populate_vertex_table() {

	// Ids of selected points
	vtkSmartPointer<vtkDataArray> ids = this -> selected_polydata -> GetPointData() -> GetArray("ids");
	for (int row = 0; row < this -> selected_polydata -> GetNumberOfPoints(); ++row) {

		// The ids of each selected point is added in the first column
		QTableWidgetItem * id_item = new QTableWidgetItem(tr("%1").arg(* (ids -> GetTuple (row))));

		table -> setItem(row, 0, id_item);

		// The coordinates of each selected point are added to the other columns
		double p[3];
		this -> selected_polydata -> GetPoint(row, p);

		QTableWidgetItem * x_item = new QTableWidgetItem(tr("%1").arg(p[0]));
		QTableWidgetItem * y_item = new QTableWidgetItem(tr("%1").arg(p[1]));
		QTableWidgetItem * z_item = new QTableWidgetItem(tr("%1").arg(p[2]));

		table -> setItem(row, 1, x_item);
		table -> setItem(row, 2, y_item);
		table -> setItem(row, 3, z_item);
	}
}

void ModifyAreaWidget::find_blob_center() {

	double blob_average_position[3]; // selected blob average coordinates
	double blob_center_position[3]; // coordinates of vertex that's closest to the blob average

	vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
	    vtkSmartPointer<vtkCenterOfMass>::New();
	centerOfMassFilter -> SetInputData(this -> selected_polydata);
	centerOfMassFilter -> SetUseScalarsAsWeights(false);
	centerOfMassFilter -> Update();
	centerOfMassFilter -> GetCenter(blob_average_position);


	// blob_average_position now contains the average location of the selected points. Note that this does not correspond to
	// a point physically present in the blob. The next step thus consists in finding the selected point that is closest to this average.
	// this point physically present in the blob is thus called the "blob center"

	unsigned int blob_center_id = this -> selected_polydata -> FindPoint(blob_average_position);
	this -> selected_polydata -> GetPoint(blob_center_id,
	                                      blob_center_position);

	// The computed locations are stored in arma::vec form
	this -> blob_center_position = {
		blob_center_position[0],
		blob_center_position[1],
		blob_center_position[2]
	};


}

void ModifyAreaWidget::update_view(int pos) {

	this -> selected_polydata -> GetPolys() -> InitTraversal();

	std::set<unsigned int> visited_vertex_indices;

	for (unsigned int selected_facet_index = 0;
	        selected_facet_index < this -> selected_polydata -> GetNumberOfPolys();
	        ++selected_facet_index) {

		vtkSmartPointer<vtkIdList> vertices_indices = vtkSmartPointer<vtkIdList>::New();

		// Current cell
		this -> selected_polydata -> GetPolys() -> GetNextCell(vertices_indices);

		// For earch vertex in the cell
		for (unsigned int local_vertex_index = 0; local_vertex_index < 3; ++ local_vertex_index) {
			unsigned int v_index = vertices_indices -> GetId(local_vertex_index);


			// If this vertex has not been previously visited
			if (visited_vertex_indices.find(v_index) == visited_vertex_indices.end()) {
				visited_vertex_indices.insert(v_index);

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

				// This point is queried
				this -> selected_polydata_original -> GetPoint(v_index, p);

				// its position is stored in a arma::vec to enable upcoming operations
				old_p_vec = {p[0], p[1], p[2]};

				switch (this -> transform_direction) {
				case TransformDirection::RADIAL :
					u = arma::normalise(old_p_vec);
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
					normalizing_constant =  this -> selected_polydata -> GetLength();
					break;

				case TransformSelection::NCLOSEST:
					normalizing_constant =  this -> active_selected_points_polydata -> GetLength();
					break;
				}


				// An interpolation factor is set depending on
				// the current choice of interpolating type.
				switch (this -> interpolation_type ) {
				case InterpolationType::UNIFORM:
					interpolating_factor = 0.1 * this -> parent -> get_asteroid()
					                       -> get_polydata() -> GetLength();
					break;

				case InterpolationType::LINEAR:
					interpolating_factor = std::max(0., 1 - arma::norm(this -> blob_center_position - old_p_vec)
					                                / normalizing_constant  ) * 0.1 * this -> parent -> get_asteroid()
					                       -> get_polydata() -> GetLength();

					break;
				case InterpolationType::PARABOLIC:
					interpolating_factor = std::max(0., 1 - std::pow(arma::norm(this -> blob_center_position - old_p_vec)
					                                / normalizing_constant , 2)) * 0.1 * this -> parent -> get_asteroid()
					                       -> get_polydata() -> GetLength();
					break;
				}

				// If the vertex is on the boundary of the selected polydata, it is clamped
				if (this -> boundary_vertex_ids_list
				        -> IsId (* this -> selected_vertices_global_ids_from_local_ids
				                 -> GetTuple(v_index) ) != -1) {
					interpolating_factor = 0;
				}

				// The position of the transform vertex is computed
				new_p_vec = old_p_vec + u * float(pos) / 100 *  interpolating_factor ;

				new_p[0] = new_p_vec(0);
				new_p[1] = new_p_vec(1);
				new_p[2] = new_p_vec(2);
				this -> selected_polydata -> GetPoints() -> SetPoint(v_index, new_p);

			}

		}

	}

	this -> selected_polydata -> Modified();
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}

void ModifyAreaWidget::accept() {

	vtkSmartPointer<vtkAppendPolyData> appendFilter =
	    vtkSmartPointer<vtkAppendPolyData>::New();

	// The normals of the selected polydata are recomputed
	vtkSmartPointer<vtkPolyDataNormals> filter = vtkPolyDataNormals::New();
	filter -> ComputePointNormalsOff();
	filter -> ComputeCellNormalsOn();
	filter -> SetInputData(this -> selected_polydata);
	filter -> Update ();
	appendFilter -> AddInputData(filter -> GetOutput());
	appendFilter -> AddInputData(this -> unselected_polydata);
	appendFilter -> Update();

	// An IDFilter is used, so as to renumber all the features
	// in a consistent way
	vtkSmartPointer<vtkIdFilter>  id_filter = vtkSmartPointer<vtkIdFilter>::New();
	id_filter -> SetIdsArrayName("ids");
	id_filter -> SetInputConnection(appendFilter -> GetOutputPort());
	id_filter -> PointIdsOn();
	id_filter -> CellIdsOn();

	id_filter -> Update();


	// The polydata is translated so as to have its coordinates
	// expressed with respect to its barycenter (constant density
	// is assumed here)

	vtkSmartPointer<vtkCenterOfMass> center_of_mass_filter =
	    vtkSmartPointer<vtkCenterOfMass>::New();

	center_of_mass_filter -> SetInputConnection(id_filter -> GetOutputPort());
	center_of_mass_filter -> SetUseScalarsAsWeights(false);
	center_of_mass_filter -> Update();

	double center[3];
	center_of_mass_filter -> GetCenter(center);

	vtkSmartPointer<vtkTransform> translation =
	    vtkSmartPointer<vtkTransform>::New();
	translation -> Translate( - center[0],  - center[1],  - center[2]);

	vtkSmartPointer<vtkTransformPolyDataFilter> translation_filter =
	    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	translation_filter -> SetInputConnection(id_filter -> GetOutputPort());
	translation_filter -> SetTransform(translation);
	translation_filter -> Update();


	this -> parent -> get_asteroid() -> get_polydata() -> DeepCopy(translation_filter -> GetOutput());
	this -> close();
}

void ModifyAreaWidget::reject() {
	this -> close();
}


void ModifyAreaWidget::show_new_slider_neighbors_pos(int pos) {
	this -> slider_neighbors_value -> setText( QString::number(pos));
}

void ModifyAreaWidget::set_new_slider_neighbors_pos() {
	this -> slider_neighbors -> setValue(this -> slider_neighbors_value -> text().toInt());

}

void ModifyAreaWidget::show_new_slider_magnitude_pos(int pos) {
	this -> slider_magnitude_value -> setText( QString::number(pos));
}

void ModifyAreaWidget::set_new_slider_magnitude_pos() {
	this -> slider_magnitude -> setValue(this -> slider_magnitude_value -> text().toInt());
}

void ModifyAreaWidget::close() {

	this -> parent -> get_renderer() -> RemoveActor(selected_actor);
	this -> parent -> get_renderer() -> RemoveActor(unselected_actor);

	this -> parent -> lateral_dockwidget -> hide();

	this -> parent -> set_actors_visibility(true);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

	QDialog::close();


}

void ModifyAreaWidget::compute_selected_cells_average_normals() {


	// The direction of the normals is averaged
	double average_normal_direction[3];
	for (int i = 0; i < this -> selected_polydata
	        -> GetCellData() -> GetNormals() -> GetNumberOfTuples(); ++i) {
		double average_normal_direction_buffer[3];

		this -> selected_polydata -> GetCellData() -> GetNormals() -> GetTuple(i,
		        average_normal_direction_buffer);

		average_normal_direction[0] = average_normal_direction[0] + average_normal_direction_buffer[0];
		average_normal_direction[1] = average_normal_direction[1] + average_normal_direction_buffer[1];
		average_normal_direction[2] = average_normal_direction[2] + average_normal_direction_buffer[2];
	}

	average_normal_direction[0] = average_normal_direction[0] / this -> selected_polydata
	                              -> GetCellData() -> GetNormals() -> GetNumberOfTuples();
	average_normal_direction[1] = average_normal_direction[1] / this -> selected_polydata
	                              -> GetCellData() -> GetNormals() -> GetNumberOfTuples();
	average_normal_direction[2] = average_normal_direction[2] / this -> selected_polydata
	                              -> GetCellData() -> GetNormals() -> GetNumberOfTuples();

	// The averaged normal is stored for later reuse
	this -> averaged_normal_array -> SetTuple(0, average_normal_direction);

}

void ModifyAreaWidget::find_N_neighbors_indices(const int N) {
	double center_point[3];
	center_point[0] = this -> blob_center_position(0);
	center_point[1] = this -> blob_center_position(1);
	center_point[2] = this -> blob_center_position(2);

	this -> point_locator -> FindClosestNPoints(N, center_point, this -> N_closest_vertices_indices);
	this -> active_selected_points_polydata -> GetPoints() -> Initialize();
	this -> selected_polydata -> GetPoints() -> GetPoints(this -> N_closest_vertices_indices, this -> active_selected_points_polydata -> GetPoints());
}

void ModifyAreaWidget::set_transform_direction(const int item_index) {
	switch (item_index) {
	case 0:
		transform_direction = TransformDirection::RADIAL;
		break;
	case 1:
		transform_direction = TransformDirection::NORMAL_AVERAGED;
		break;

	default:
		std::cout << " Case not implemented in set_transform_direction. Got item_index== " << item_index << std::endl;
		break;

	}
}

void ModifyAreaWidget::set_interpolation_type(const int item_index) {
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

void ModifyAreaWidget::set_transform_selection(const int item_index) {
	switch (item_index) {
	case 0:
		transform_selection = TransformSelection::SELECTED;
		this -> slider_neighbors_title -> setVisible(false);
		this -> slider_neighbors_holder_widget -> setVisible(false);
		this -> slider_neighbors -> setValue(this -> selected_polydata -> GetNumberOfPoints());
		// The line above is necessary because it will call the slot updating the view

		break;
	case 1:
		transform_selection = TransformSelection::NCLOSEST;
		this -> slider_neighbors_title -> setVisible(true);
		this -> slider_neighbors_holder_widget -> setVisible(true);
		this -> slider_neighbors -> setMaximum(this -> selected_polydata -> GetNumberOfPoints());
		this -> slider_neighbors -> setValue(this -> selected_polydata -> GetNumberOfPoints());
		break;
	default:
		std::cout << " Case not implemented in set_transform_selection. Got item_index== " << item_index << std::endl;
		break;

	}
}

void ModifyAreaWidget::update_view_unchanged_slider_magnitude() {
	int pos = this -> slider_magnitude -> value();
	this -> update_view(pos);
}

void ModifyAreaWidget::find_N_neighbors_indices_and_update_view(const int N) {
	this -> find_N_neighbors_indices(N);
	this -> update_view_unchanged_slider_magnitude();
}