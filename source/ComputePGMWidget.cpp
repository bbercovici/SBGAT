#include "ComputePGMWidget.hpp"

ComputePGMWidget::ComputePGMWidget(Mainwindow * parent) {
	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;
	this -> parent -> set_action_status(false, this -> parent -> openComputePGMWidgetAct);

	this -> main_layout = new QVBoxLayout();
	this -> physical_properties_layout = new QGridLayout();
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, this);
	this ->	compute_PGM_button = new QPushButton("Compute Polyhedron Gravity Model", this);
	this -> physical_properties_box = new QGroupBox("Physical Properties", this);

	this -> density_title_label = new QLabel("Density: ", this);
	this -> density_unit_label = new QLabel("kg/m^3", this);
	this -> density_qlineedit = new QLineEdit(this);
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



	this -> main_layout -> addWidget(this -> physical_properties_box);
	this -> main_layout -> addWidget(this -> compute_PGM_button);
	this -> main_layout -> addWidget(this -> button_box);

	this -> main_layout -> addStretch(1);
	this -> setLayout(this -> main_layout);


	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(close()));
	connect(this -> compute_PGM_button, SIGNAL(clicked()),
	        this, SLOT(compute_pgm()));
}

void ComputePGMWidget::close() {


	this -> parent -> lateral_dockwidget -> hide();
	this -> parent -> set_action_status(true, this -> parent -> openComputePGMWidgetAct);


	QDialog::close();
}

void ComputePGMWidget::compute_pgm() {

	// InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	//         -> GetInteractorStyle());
	// vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();

	// Asteroid asteroid = Asteroid(all_points_polydata,)

}
