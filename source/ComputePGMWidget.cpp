#include "ComputePGMWidget.hpp"

ComputePGMWidget::ComputePGMWidget(Mainwindow * parent) {
	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;

	this -> main_layout = new QVBoxLayout();
	this -> physical_properties_layout = new QGridLayout();
	this -> compute_acceleration_layout = new QGridLayout();
	this -> compute_pgm_layout = new QGridLayout();
	this -> visualize_pgm_layout = new QGridLayout();

	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, this);
	this ->	compute_local_acceleration_button = new QPushButton("Compute Local Acceleration", this);
	this ->	compute_global_acceleration_button = new QPushButton("Compute Surface Acceleration", this);

	this -> show_acceleration_magnitude_button = new QPushButton("Acceleration magnitude", this);
	this -> show_normal_acceleration_angle_button = new QPushButton("Normal/acceleration angle", this);
	this -> show_radial_acceleration_angle_button = new QPushButton("Radial/acceleration angle", this);
	this -> show_radial_acceleration_component_button = new QPushButton("Radial acceleration component", this);
	this -> show_normal_acceleration_component_button = new QPushButton("Normal acceleration component", this);;
	this -> show_orthonormal_acceleration_magnitude_button = new QPushButton("Orthonormal acceleration magnitude", this);
	this -> show_orthoradial_acceleration_magnitude_button = new QPushButton("Orthoradial acceleration magnitude", this);



	this -> physical_properties_box = new QGroupBox("Physical Properties", this);
	this ->	compute_acceleration_box = new QGroupBox("Local Gravity Acceleration", this);
	this -> compute_pgm_box = new QGroupBox("Surface Gravity Acceleration");
	this -> visualize_pgm_box = new QGroupBox("Visualization");

	this -> density_title_label = new QLabel("Density: ", this);
	this -> density_unit_label = new QLabel("kg/m^3", this);

	this -> x_coordinate_label = new QLabel("x (m) ", this);
	this -> y_coordinate_label = new QLabel("y (m) ", this);
	this -> z_coordinate_label = new QLabel("z (m) ", this);

	this -> x_coordinate_qlineedit = new QLineEdit(this);
	this -> y_coordinate_qlineedit = new QLineEdit(this);
	this -> z_coordinate_qlineedit = new QLineEdit(this);

	this -> x_coordinate_qlineedit -> setValidator( new QDoubleValidator(-100, 100, 10, this) );
	this -> y_coordinate_qlineedit -> setValidator( new QDoubleValidator(-100, 100, 10, this) );
	this -> z_coordinate_qlineedit -> setValidator( new QDoubleValidator(-100, 100, 10, this) );

	this -> compute_acceleration_plainedit = new QPlainTextEdit(this);
	this -> compute_acceleration_plainedit -> setReadOnly(true);

	this -> density_qlineedit = new QLineEdit(this);
	this -> density_qlineedit -> setText(QString::number(1e3));
	this -> density_qlineedit -> setValidator( new QDoubleValidator(0, 100, 10, this) );

	this -> gravitational_constant_title_label = new QLabel("Gravitational Constant: ", this);
	this -> gravitational_constant_unit_label = new QLabel("m^3/ (kg s^2)", this);
	this -> gravitational_constant_qlineedit = new QLineEdit(this);
	this -> gravitational_constant_qlineedit -> setText(QString::number(arma::datum::G));
	this -> gravitational_constant_qlineedit -> setValidator( new QDoubleValidator(0, 100, 10, this) );

	this -> physical_properties_box -> setLayout(this -> physical_properties_layout);
	this -> physical_properties_layout -> addWidget(this -> density_title_label, 0, 0, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> density_qlineedit, 0, 1, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> density_unit_label, 0, 2, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_title_label, 1, 0, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_qlineedit, 1, 1, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_unit_label, 1, 2, 1, 1);

	this -> compute_acceleration_box -> setLayout(this -> compute_acceleration_layout);
	this -> compute_acceleration_layout -> addWidget(this -> x_coordinate_label, 0, 0, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> y_coordinate_label, 0, 1, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> z_coordinate_label, 0, 2, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> x_coordinate_qlineedit, 1, 0, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> y_coordinate_qlineedit, 1, 1, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> z_coordinate_qlineedit, 1, 2, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> compute_local_acceleration_button, 2, 0, 1, 3);
	this -> compute_acceleration_layout -> addWidget(this -> compute_acceleration_plainedit, 3, 0, 2, 3);


	this -> compute_pgm_box -> setLayout(this -> compute_pgm_layout);
	this -> compute_pgm_layout -> addWidget(this -> compute_global_acceleration_button, 0, 0, 1, 1);
	this -> compute_pgm_layout -> addWidget(this -> visualize_pgm_box, 1, 0, 1, 1);

	this -> visualize_pgm_box -> setLayout(this -> visualize_pgm_layout);
	this -> visualize_pgm_layout -> addWidget(this -> show_acceleration_magnitude_button, 0, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_normal_acceleration_angle_button, 1, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_radial_acceleration_angle_button, 2, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_radial_acceleration_component_button, 3, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_normal_acceleration_component_button, 4, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_orthoradial_acceleration_magnitude_button, 5, 0, 1, 1);
	this -> visualize_pgm_layout -> addWidget(this -> show_orthonormal_acceleration_magnitude_button, 6, 0, 1, 1);

	this -> visualize_pgm_box -> setDisabled(1);

	this -> main_layout -> addWidget(this -> physical_properties_box);
	this -> main_layout -> addWidget(this -> compute_acceleration_box);
	this -> main_layout -> addWidget(this -> compute_pgm_box);

	this -> main_layout -> addWidget(this -> button_box);

	this -> main_layout -> addStretch(1);
	this -> setLayout(this -> main_layout);

	// Creates a sphere to represent specific location
	// at which the gravity field must be evaluated

	InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	        -> GetInteractorStyle());
	vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();

	// A OBBTree is constructed in order to check
	// if a specified point is inside or outside the shape


	this -> shape_mod_obb_tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
	this -> shape_mod_obb_tree -> SetDataSet(all_points_polydata);
	this -> shape_mod_obb_tree -> BuildLocator();


	this -> specific_point_sphere =
	    vtkSmartPointer<vtkSphereSource>::New();

	this -> specific_point_sphere -> SetCenter(0.0, 0.0, 0.0);
	this -> specific_point_sphere -> SetPhiResolution(100);
	this -> specific_point_sphere -> SetThetaResolution(100);
	this -> specific_point_sphere -> SetRadius(all_points_polydata -> GetLength() / 750);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
	    vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper -> SetInputConnection( this -> specific_point_sphere -> GetOutputPort());

	this -> specific_point_actor =
	    vtkSmartPointer<vtkActor>::New();
	this -> specific_point_actor -> SetMapper(mapper);
	this -> specific_point_actor -> GetProperty() -> SetColor(0.3, 0.8, 0.2); //(R,G,B)

	this -> parent -> get_renderer() -> AddActor(this -> specific_point_actor);

	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(close()));
	connect(this -> compute_local_acceleration_button, SIGNAL(clicked()),
	        this, SLOT(compute_local_pgm()));

	connect(this -> compute_global_acceleration_button, SIGNAL(clicked()),
	        this, SLOT(compute_global_pgm()));

	connect(this -> x_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_local_acceleration_sphere()));
	connect(this -> y_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_local_acceleration_sphere()));
	connect(this -> z_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_local_acceleration_sphere()));

	connect(this -> show_acceleration_magnitude_button, SIGNAL(clicked()),
	        this, SLOT(show_acceleration_magnitude()));

	connect(this -> show_normal_acceleration_angle_button, SIGNAL(clicked()),
	        this, SLOT(show_normal_acceleration_angle()));

	connect(this -> show_radial_acceleration_angle_button, SIGNAL(clicked()),
	        this, SLOT(show_radial_acceleration_angle()));

	connect(this -> show_radial_acceleration_component_button, SIGNAL(clicked()),
	        this, SLOT(show_radial_acceleration_component()));

	connect(this -> show_normal_acceleration_component_button, SIGNAL(clicked()),
	        this, SLOT(show_normal_acceleration_component()));

	connect(this -> show_orthoradial_acceleration_magnitude_button, SIGNAL(clicked()),
	        this, SLOT(show_orthoradial_acceleration_magnitude()));

	connect(this -> show_orthonormal_acceleration_magnitude_button, SIGNAL(clicked()),
	        this, SLOT(show_orthonormal_acceleration_magnitude()));



	this -> parent -> statusBar() -> showMessage("Constructing PGM");

	double G_dens = (this -> density_qlineedit -> text().toDouble())
	                * (this -> gravitational_constant_qlineedit -> text().toDouble());
	this -> asteroid = new Asteroid(all_points_polydata, G_dens);
	this -> parent -> statusBar() -> showMessage("Ready");


}

void ComputePGMWidget::move_local_acceleration_sphere() {
	double x = this -> x_coordinate_qlineedit -> text().toDouble();
	double y = this -> y_coordinate_qlineedit -> text().toDouble();
	double z = this -> z_coordinate_qlineedit -> text().toDouble();

	this -> specific_point_sphere -> SetCenter(x, y, z);
	this -> specific_point_sphere -> Update();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();



}

void ComputePGMWidget::close() {
	this -> parent -> get_renderer() -> RemoveActor(this -> specific_point_actor);

	this -> parent -> lateral_dockwidget -> hide();
	this -> parent -> set_action_status(true, this -> parent -> openComputePGMWidgetAct);
	delete(asteroid);

	this -> cleanup();

	QDialog::close();
}

void ComputePGMWidget::compute_local_pgm() {
	this -> parent -> statusBar() -> showMessage("Computing acceleration");

	this -> compute_acceleration_plainedit -> clear();

	InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	        -> GetInteractorStyle());
	vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();

	// A ray is cast in the outward direction to check if the specified
	// point is inside the polyhedron
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();


	double specified_point[3] = {
		x_coordinate_qlineedit -> text().toDouble(),
		y_coordinate_qlineedit -> text().toDouble(),
		z_coordinate_qlineedit -> text().toDouble()
	};

	arma::vec specified_point_arma = {
		specified_point[0],
		specified_point[1],
		specified_point[2]
	};

	arma::vec arbitrary_end_point_arma = specified_point_arma * (1 + 2 * all_points_polydata -> GetLength() / (arma::norm(specified_point_arma)));

	double arbitrary_end_point[3] = {
		arbitrary_end_point_arma(0),
		arbitrary_end_point_arma(1),
		arbitrary_end_point_arma(2)

	};

	int hit = this -> shape_mod_obb_tree -> IntersectWithLine(
	              specified_point,
	              arbitrary_end_point,
	              1e-6,
	              points,
	              cellIds);


	if (hit != 0 || arma::norm(specified_point_arma) < 1e-9) {
		this -> compute_acceleration_plainedit
		-> appendPlainText(
		    QString::fromStdString("The specified point must be outside of the polyhedron"));
	}

	else {

		Vect Xsc(3);

		Xsc[0] = specified_point[0];
		Xsc[1] = specified_point[1];
		Xsc[2] = specified_point[2];


		double G_dens = (this -> density_qlineedit -> text().toDouble())
		                * (this -> gravitational_constant_qlineedit -> text().toDouble());
		this -> asteroid -> setmGs(G_dens);

		Vect grav = this -> asteroid -> PolyGrav(Xsc, true);

		this -> compute_acceleration_plainedit -> appendPlainText(
		    QString::fromStdString("Acceleration at specified point (m/s^2): "));
		this -> compute_acceleration_plainedit -> appendPlainText("x:  " + QString::number(grav[0]));
		this -> compute_acceleration_plainedit -> appendPlainText("y:  " + QString::number(grav[1]));
		this -> compute_acceleration_plainedit -> appendPlainText("z:  " + QString::number(grav[2]));
	}
	this -> parent -> statusBar() -> showMessage("Ready");


}

void ComputePGMWidget::compute_global_pgm() {

	double G_dens = (this -> density_qlineedit -> text().toDouble())
	                * (this -> gravitational_constant_qlineedit -> text().toDouble());
	this -> asteroid -> setmGs(G_dens);
	this -> asteroid -> compute_global_pgm();
	this -> visualize_pgm_box -> setEnabled(1);

}

void ComputePGMWidget::show_acceleration_magnitude() {
	this -> cleanup();

	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_mag = - std::numeric_limits<double>::infinity();
	double min_mag = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};

		if (arma::norm(acceleration) > max_mag) {
			max_mag = arma::norm(acceleration) ;
		}

		if (arma::norm(acceleration) < min_mag) {
			min_mag = arma::norm(acceleration) ;
		}

		surface_data -> InsertNextValue(arma::norm(acceleration));
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_mag, max_mag);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Acceleration Magnitude (m/s^2)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

void ComputePGMWidget::cleanup() {
	if (this -> parent -> get_renderer() -> GetActors2D() -> GetNumberOfItems() != 0) {
		this -> parent -> get_renderer() -> RemoveActor2D(
		    this -> parent -> get_renderer() -> GetActors2D() -> GetLastActor2D());
	}
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];
	shape_actor -> GetMapper() -> ScalarVisibilityOff();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}

void ComputePGMWidget::show_normal_acceleration_angle() {
	this -> cleanup();

	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_angle = - std::numeric_limits<double>::infinity();
	double min_angle = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};

		arma::vec normal = {
			this -> asteroid -> get_ListN()[i][0],
			this -> asteroid -> get_ListN()[i][1],
			this -> asteroid -> get_ListN()[i][2]
		};

		double angle = 180. / arma::datum::pi * std::acos(arma::dot(acceleration, normal) / arma::norm(acceleration));

		if (angle > max_angle) {
			max_angle = angle ;
		}

		if (angle < min_angle) {
			min_angle = angle ;
		}

		surface_data -> InsertNextValue(angle);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_angle, max_angle);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Normal/acceleration angle (deg)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

void ComputePGMWidget::show_radial_acceleration_angle() {
	this -> cleanup();

	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_angle = - std::numeric_limits<double>::infinity();
	double min_angle = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};
		unsigned int P1_index = this -> asteroid -> get_ListTri()[i][0];
		unsigned int P2_index = this -> asteroid -> get_ListTri()[i][1];
		unsigned int P3_index = this -> asteroid -> get_ListTri()[i][2];

		arma::vec P1 = {
			this -> asteroid -> get_X()[P1_index],
			this -> asteroid -> get_Y()[P1_index],
			this -> asteroid -> get_Z()[P1_index]
		};

		arma::vec P2 = {
			this -> asteroid -> get_X()[P2_index],
			this -> asteroid -> get_Y()[P2_index],
			this -> asteroid -> get_Z()[P2_index]
		};

		arma::vec P3 = {
			this -> asteroid -> get_X()[P3_index],
			this -> asteroid -> get_Y()[P3_index],
			this -> asteroid -> get_Z()[P3_index]
		};

		arma::vec P = 1. / 3. * (P1 + P2 + P3);
		arma::vec radial = P / arma::norm(P);

		double angle = 180. / arma::datum::pi * std::acos(arma::dot(acceleration, radial) / arma::norm(acceleration));

		if (angle > max_angle) {
			max_angle = angle ;
		}

		if (angle < min_angle) {
			min_angle = angle ;
		}

		surface_data -> InsertNextValue(angle);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_angle, max_angle);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Radial/acceleration angle (deg)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}

void ComputePGMWidget::show_radial_acceleration_component() {
	this -> cleanup();


	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_component = - std::numeric_limits<double>::infinity();
	double min_component = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};
		unsigned int P1_index = this -> asteroid -> get_ListTri()[i][0];
		unsigned int P2_index = this -> asteroid -> get_ListTri()[i][1];
		unsigned int P3_index = this -> asteroid -> get_ListTri()[i][2];

		arma::vec P1 = {
			this -> asteroid -> get_X()[P1_index],
			this -> asteroid -> get_Y()[P1_index],
			this -> asteroid -> get_Z()[P1_index]
		};

		arma::vec P2 = {
			this -> asteroid -> get_X()[P2_index],
			this -> asteroid -> get_Y()[P2_index],
			this -> asteroid -> get_Z()[P2_index]
		};

		arma::vec P3 = {
			this -> asteroid -> get_X()[P3_index],
			this -> asteroid -> get_Y()[P3_index],
			this -> asteroid -> get_Z()[P3_index]
		};
		this -> asteroid -> get_X()[P3_index] = 1;

		arma::vec P = 1. / 3. * (P1 + P2 + P3);
		arma::vec radial = P / arma::norm(P);

		double component = arma::dot(acceleration, radial) ;

		if (component > max_component) {
			max_component = component ;
		}

		if (component < min_component) {
			min_component = component ;
		}

		surface_data -> InsertNextValue(component);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_component, max_component);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();
	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Radial acceleration component (m/s^2)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();



}

void ComputePGMWidget::show_normal_acceleration_component() {
	this -> cleanup();


	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_component = - std::numeric_limits<double>::infinity();
	double min_component = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};

		arma::vec normal = {
			this -> asteroid -> get_ListN()[i][0],
			this -> asteroid -> get_ListN()[i][1],
			this -> asteroid -> get_ListN()[i][2]
		};

		double component = arma::dot(acceleration, normal) ;

		if (component > max_component) {
			max_component = component ;
		}

		if (component < min_component) {
			min_component = component ;
		}

		surface_data -> InsertNextValue(component);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_component, max_component);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();

	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Normal acceleration component (m/s^2)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

void ComputePGMWidget::show_orthoradial_acceleration_magnitude() {
	this -> cleanup();

	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_magnitude = - std::numeric_limits<double>::infinity();
	double min_magnitude = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};

		unsigned int P1_index = this -> asteroid -> get_ListTri()[i][0];
		unsigned int P2_index = this -> asteroid -> get_ListTri()[i][1];
		unsigned int P3_index = this -> asteroid -> get_ListTri()[i][2];

		arma::vec P1 = {
			this -> asteroid -> get_X()[P1_index],
			this -> asteroid -> get_Y()[P1_index],
			this -> asteroid -> get_Z()[P1_index]
		};

		arma::vec P2 = {
			this -> asteroid -> get_X()[P2_index],
			this -> asteroid -> get_Y()[P2_index],
			this -> asteroid -> get_Z()[P2_index]
		};

		arma::vec P3 = {
			this -> asteroid -> get_X()[P3_index],
			this -> asteroid -> get_Y()[P3_index],
			this -> asteroid -> get_Z()[P3_index]
		};

		arma::vec P = 1. / 3. * (P1 + P2 + P3);
		arma::vec radial = P / arma::norm(P);

		double magnitude = arma::norm(acceleration
		                              - arma::dot(acceleration, radial) * radial) ;

		if (magnitude > max_magnitude) {
			max_magnitude = magnitude ;
		}

		if (magnitude < min_magnitude) {
			min_magnitude = magnitude ;
		}

		surface_data -> InsertNextValue(magnitude);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_magnitude, max_magnitude);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();

	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Orthoradial acceleration magnitude (m/s^2)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();




}

void ComputePGMWidget::show_orthonormal_acceleration_magnitude() {
	this -> cleanup();


	// Create cell data
	vtkSmartPointer<vtkDoubleArray> surface_data =
	    vtkSmartPointer<vtkDoubleArray>::New();


	double max_magnitude = - std::numeric_limits<double>::infinity();
	double min_magnitude = std::numeric_limits<double>::infinity();
	for (unsigned int i = 0; i < this -> asteroid -> GetNOF(); i++) {

		arma::vec acceleration = {
			this -> asteroid -> get_surface_grav()[i][0],
			this -> asteroid -> get_surface_grav()[i][1],
			this -> asteroid -> get_surface_grav()[i][2]
		};

		arma::vec normal = {
			this -> asteroid -> get_ListN()[i][0],
			this -> asteroid -> get_ListN()[i][1],
			this -> asteroid -> get_ListN()[i][2]
		};

		double magnitude = arma::norm(acceleration
		                              - arma::dot(acceleration, normal) * normal) ;

		if (magnitude > max_magnitude) {
			max_magnitude = magnitude ;
		}

		if (magnitude < min_magnitude) {
			min_magnitude = magnitude ;
		}

		surface_data -> InsertNextValue(magnitude);
	}

	// Create a lookup table to map cell data to colors
	vtkSmartPointer<vtkLookupTable> lut =
	    vtkSmartPointer<vtkLookupTable>::New();

	lut -> SetNumberOfTableValues(this -> asteroid -> GetNOF());
	lut -> Build();

	//
	vtkActor * shape_actor = this -> parent -> get_actor_vector()[0];

	shape_actor -> GetMapper() -> SetScalarRange(min_magnitude, max_magnitude);

	shape_actor -> GetMapper() -> SetLookupTable(lut);
	shape_actor -> GetMapper() -> GetInput() -> GetCellData() -> SetScalars(surface_data);
	shape_actor -> GetMapper() -> ScalarVisibilityOn();
	shape_actor -> GetMapper() -> SetScalarModeToUseCellData();

	vtkSmartPointer<vtkScalarBarActor> scalarBar =
	    vtkSmartPointer<vtkScalarBarActor>::New();

	scalarBar -> SetLookupTable(shape_actor -> GetMapper() -> GetLookupTable());
	scalarBar -> SetTitle("Orthonormal acceleration magnitude (m/s^2)");
	scalarBar -> SetUnconstrainedFontSize (true);
	scalarBar -> GetTitleTextProperty() -> SetFontSize(10);
	scalarBar -> GetLabelTextProperty() -> SetFontSize(10);
	scalarBar -> SetNumberOfLabels(4);


	this -> parent -> get_renderer() -> AddActor2D(scalarBar);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();



}

