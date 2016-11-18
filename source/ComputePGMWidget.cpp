#include "ComputePGMWidget.hpp"

ComputePGMWidget::ComputePGMWidget(Mainwindow * parent) {
	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;

	this -> main_layout = new QVBoxLayout();
	this -> physical_properties_layout = new QGridLayout();
	this -> compute_acceleration_layout = new QGridLayout();

	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, this);
	this ->	compute_PGM_button = new QPushButton("Compute Polyhedron Gravity Model", this);
	this -> physical_properties_box = new QGroupBox("Physical Properties", this);
	this ->	compute_acceleration_box = new QGroupBox("Gravity acceleration", this);

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

	this -> scaling_factor_title_label = new QLabel("Scaling Factor: ", this);
	this -> scaling_factor_qlineedit = new QLineEdit(this);
	this -> scaling_factor_qlineedit -> setText(QString::number(1));
	this -> scaling_factor_qlineedit -> setValidator( new QDoubleValidator(0, 100, 10, this) );

	this -> physical_properties_box -> setLayout(this -> physical_properties_layout);
	this -> physical_properties_layout -> addWidget(this -> density_title_label, 0, 0, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> density_qlineedit, 0, 1, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> density_unit_label, 0, 2, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_title_label, 1, 0, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_qlineedit, 1, 1, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> gravitational_constant_unit_label, 1, 2, 1, 1);

	this -> physical_properties_layout -> addWidget(this -> scaling_factor_title_label, 2, 0, 1, 1);
	this -> physical_properties_layout -> addWidget(this -> scaling_factor_qlineedit, 2, 1, 1, 1);

	this -> compute_acceleration_box -> setLayout(this -> compute_acceleration_layout);
	this -> compute_acceleration_layout -> addWidget(this -> x_coordinate_label, 0, 0, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> y_coordinate_label, 0, 1, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> z_coordinate_label, 0, 2, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> x_coordinate_qlineedit, 1, 0, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> y_coordinate_qlineedit, 1, 1, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> z_coordinate_qlineedit, 1, 2, 1, 1);
	this -> compute_acceleration_layout -> addWidget(this -> compute_PGM_button, 2, 0, 1, 3);
	this -> compute_acceleration_layout -> addWidget(this -> compute_acceleration_plainedit, 3, 0, 2, 3);

	this -> main_layout -> addWidget(this -> physical_properties_box);
	this -> main_layout -> addWidget(this -> compute_acceleration_box);

	this -> main_layout -> addWidget(this -> button_box);

	this -> main_layout -> addStretch(1);
	this -> setLayout(this -> main_layout);

	// Creates a sphere to represent specific location
	// at which the gravity field must be evaluated

	InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	        -> GetInteractorStyle());
	vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();


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
	connect(this -> compute_PGM_button, SIGNAL(clicked()),
	        this, SLOT(compute_pgm()));
	connect(this -> x_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_specific_point_sphere()));
	connect(this -> y_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_specific_point_sphere()));
	connect(this -> z_coordinate_qlineedit, SIGNAL(textChanged(QString)),
	        this, SLOT(move_specific_point_sphere()));
}

void ComputePGMWidget::move_specific_point_sphere() {
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


	QDialog::close();
}

void ComputePGMWidget::compute_pgm() {

	this -> compute_acceleration_plainedit -> clear();

	InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	        -> GetInteractorStyle());
	vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();

	// A ray is cast in the outward direction to check if the specified
	// point is inside the polyhedron
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();


	vtkSmartPointer<vtkModifiedBSPTree> shape_mod_obb_tree = vtkSmartPointer<vtkModifiedBSPTree>::New();
	shape_mod_obb_tree -> SetDataSet(all_points_polydata);
	shape_mod_obb_tree -> BuildLocator();

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

	int hit = shape_mod_obb_tree -> IntersectWithLine(
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

		double G_dens = (this -> density_qlineedit -> text().toDouble())
		                * (this -> gravitational_constant_qlineedit -> text().toDouble());
		Asteroid asteroid = Asteroid(all_points_polydata, G_dens);

		Vect Xsc(3);

		Xsc[0] = specified_point[0];
		Xsc[1] = specified_point[1];
		Xsc[2] = specified_point[2];

		Vect grav = asteroid.PolyGrav(Xsc, true);

		this -> compute_acceleration_plainedit -> appendPlainText(
		    QString::fromStdString("Acceleration at specified point (m/s^2): "));
		this -> compute_acceleration_plainedit -> appendPlainText("x:  " + QString::number(grav[0]));
		this -> compute_acceleration_plainedit -> appendPlainText("y:  " + QString::number(grav[1]));
		this -> compute_acceleration_plainedit -> appendPlainText("z:  " + QString::number(grav[2]));
	}

}
