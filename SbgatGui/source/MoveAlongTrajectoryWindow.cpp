#include "MoveAlongTrajectoryWindow.hpp"
#include <QVBoxLayout>
#include <QDialogButtonBox>

#include <QLabel>
#include <QVTKOpenGLWidget.h>
#include <vtkCamera.h>

#include <armadillo>


using namespace SBGAT_GUI;

MoveAlongTrajectoryWindow::MoveAlongTrajectoryWindow(Mainwindow * parent) : QDialog(parent,Qt::WindowStaysOnTopHint) {

	this -> parent = parent;
	this -> setWindowTitle("Move spacecraft along trajectory");

	QVBoxLayout * settings_layout = new QVBoxLayout(this);


	QGroupBox * trajectory_group = new QGroupBox(tr("Trajectory"));
	QGroupBox * spacecraft_group = new QGroupBox(tr("Spacecraft"));
	QGroupBox * camera_group = new QGroupBox(tr("Camera behavior"));
	this -> position_group = new QGroupBox(tr("Set position"));

	
	QGridLayout * trajectory_group_layout = new QGridLayout(trajectory_group);
	QGridLayout * spacecraft_group_layout = new QGridLayout(spacecraft_group);
	QGridLayout * camera_group_layout = new QGridLayout(camera_group);
	QGridLayout * position_group_layout = new QGridLayout(position_group);

	QLabel * trajectory_label = new QLabel("Select trajectory: ", this);
	QLabel * spacecraft_label = new QLabel("Select spacecraft: ", this);
	QLabel * camera_focus_label = new QLabel("Focus on satellite: ", this);
	QLabel * camera_orientation_label = new QLabel("Orientation: ", this);


	this -> camera_focus_check = new QCheckBox(this);


	this -> pos_spinbox = new QSpinBox(this);
	this -> trajectory_combo_box = new QComboBox (this);
	this -> spacecraft_combo_box = new QComboBox (this);
	this -> camera_orientation_combo_box = new QComboBox (this);
	this -> pos_slider = new QSlider(Qt::Horizontal,this);
	this -> pos_slider -> setRange(0, 100);

	// Spacecraft 
	spacecraft_group_layout -> addWidget(spacecraft_label, 0, 0, 1, 1);
	spacecraft_group_layout -> addWidget(this -> spacecraft_combo_box, 0, 1, 1, 1);

	trajectory_group_layout -> addWidget(trajectory_label, 0, 0, 1, 1);
	trajectory_group_layout -> addWidget(this -> trajectory_combo_box, 0, 1, 1, 1);

	camera_group_layout -> addWidget(camera_focus_label, 0, 0, 1, 1);
	camera_group_layout -> addWidget(this -> camera_focus_check, 0, 1, 1, 1);

	camera_group_layout -> addWidget(camera_orientation_label, 1, 0, 1, 1);
	camera_group_layout -> addWidget(this -> camera_orientation_combo_box, 1, 1, 1, 1);


	position_group_layout -> addWidget(this -> pos_slider, 0, 0, 1, 1);
	position_group_layout -> addWidget(this -> pos_spinbox, 0, 1, 1, 1);


	// Creating the button box
	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);

	// Connecting buttons signals to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));

	connect(this -> pos_slider,SIGNAL(valueChanged(int)),this ,SLOT(update_position()));
	connect(this -> pos_spinbox,SIGNAL(valueChanged(int)),this ,SLOT(update_position()));
	connect(this -> trajectory_combo_box,SIGNAL(currentIndexChanged(int)),this,SLOT(changed_trajectory()));
	
	connect(this -> camera_focus_check,SIGNAL(stateChanged(int)),this,SLOT(toggle_camera_focus(int)));




	settings_layout -> addWidget(trajectory_group);
	settings_layout -> addWidget(spacecraft_group);
	settings_layout -> addWidget(camera_group);
	settings_layout -> addWidget(position_group);


	this -> init();

	
}

void MoveAlongTrajectoryWindow::toggle_camera_focus(int state){

	if (state == Qt::Checked){
		this -> camera_orientation_combo_box -> setEnabled(true);
	}
	else{
		this -> camera_orientation_combo_box -> setEnabled(false);
	}

}



void MoveAlongTrajectoryWindow::init(){
	auto wrapped_trajectory_data = this -> parent -> get_wrapped_trajectory_data();
	
	for (auto it = wrapped_trajectory_data.begin(); it != wrapped_trajectory_data.end(); ++it){
		this -> trajectory_combo_box -> insertItem(this -> trajectory_combo_box -> count(),
			QString::fromStdString(it -> first));
	}

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		this -> spacecraft_combo_box -> insertItem(this -> spacecraft_combo_box -> count(),QString::fromStdString(it -> first));
		
	}

	this -> camera_orientation_combo_box -> insertItem(0,"Radial");
	this -> camera_orientation_combo_box -> insertItem(1,"Trailing");
	this -> camera_orientation_combo_box -> insertItem(1,"Orbit");


	this -> toggle_camera_focus(Qt::Unchecked);

}

void MoveAlongTrajectoryWindow::changed_trajectory(){
	
	if (this -> trajectory_combo_box -> count() > 0){
	// the slider range is adjusted to reflect the current trajectory
		std::string trajectory_name = this -> trajectory_combo_box -> currentText().toStdString();
		auto wrapped_trajectory_data = this -> parent -> get_wrapped_trajectory_data();
		unsigned int N_points = wrapped_trajectory_data[trajectory_name] -> get_points() -> GetNumberOfPoints();
		this -> pos_slider -> setRange(0, N_points - 1);
		this -> pos_spinbox -> setRange(0, N_points - 1);
	}

}

void MoveAlongTrajectoryWindow::update_position(){

	std::string spacecraft_name = this -> spacecraft_combo_box -> currentText().toStdString();
	std::string trajectory_name = this -> trajectory_combo_box -> currentText().toStdString();

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_trajectory_data = this -> parent -> get_wrapped_trajectory_data();

	unsigned int N_points = wrapped_trajectory_data[trajectory_name] -> get_points() -> GetNumberOfPoints();

	QObject* obj = sender();
	if( obj == this -> pos_slider ){ 
		this -> pos_spinbox -> setValue(this -> pos_slider -> value());
	}
	else{
		this -> pos_slider -> setValue(this -> pos_spinbox -> value());
	}

	
	auto spacecraft_actor = wrapped_spacecraft_data[spacecraft_name] -> get_actor();
	double pos[3];

	wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value(),pos);

	spacecraft_actor -> SetPosition(pos[0], pos[1], pos[2]);  
	spacecraft_actor -> Modified();



	if (this -> camera_focus_check -> checkState() == Qt::Checked){

		auto prop_to_focus_on = wrapped_spacecraft_data[spacecraft_name] -> get_actor();

		vtkSmartPointer<vtkCamera> camera = this -> parent -> get_renderer() -> GetActiveCamera();
		camera -> SetFocalPoint(prop_to_focus_on -> GetPosition());

		if (this -> camera_orientation_combo_box -> currentText() == "Radial"){
			double prop_pos[3];
			double prop_bounds[6];

			prop_to_focus_on -> GetPosition(prop_pos);
			prop_to_focus_on -> GetBounds (prop_bounds);

			arma::vec prop_pos_arma = {prop_pos[0],prop_pos[1],prop_pos[2]};
			arma::vec bbox_min = {prop_bounds[0],prop_bounds[2],prop_bounds[4]};
			arma::vec bbox_max = {prop_bounds[1],prop_bounds[3],prop_bounds[5]};
			double diagonal = arma::norm(bbox_max - bbox_min);

			arma::vec cam_pos = prop_pos_arma + 4 * diagonal * arma::normalise(prop_pos_arma);
			camera -> SetPosition(cam_pos.colptr(0));

		}
		else if (this -> camera_orientation_combo_box -> currentText() == "Trailing"){
			double prop_pos[3];
			double prop_bounds[6];

			prop_to_focus_on -> GetPosition(prop_pos);
			prop_to_focus_on -> GetBounds (prop_bounds);

			arma::vec prop_pos_arma = {prop_pos[0],prop_pos[1],prop_pos[2]};
			arma::vec bbox_min = {prop_bounds[0],prop_bounds[2],prop_bounds[4]};
			arma::vec bbox_max = {prop_bounds[1],prop_bounds[3],prop_bounds[5]};
			double diagonal = arma::norm(bbox_max - bbox_min);


			arma::vec pos_before(3);
			arma::vec pos_now(3);

			if (this -> pos_slider -> value() != 0){
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value() - 1,pos_before.colptr(0));
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value(),pos_now.colptr(0));
			}
			else{
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value(),pos_before.colptr(0));
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value() + 1,pos_now.colptr(0));

			}

			arma::vec cam_pos = prop_pos_arma + 4 * diagonal * arma::normalise(pos_before - pos_now);
			camera -> SetPosition(cam_pos.colptr(0));

		}

		else if (this -> camera_orientation_combo_box -> currentText() == "Orbit"){
			double prop_pos[3];
			double prop_bounds[6];
			arma::vec pos_before(3);
			arma::vec pos_now(3);

			prop_to_focus_on -> GetPosition(prop_pos);
			prop_to_focus_on -> GetBounds (prop_bounds);

			arma::vec prop_pos_arma = {prop_pos[0],prop_pos[1],prop_pos[2]};
			arma::vec bbox_min = {prop_bounds[0],prop_bounds[2],prop_bounds[4]};
			arma::vec bbox_max = {prop_bounds[1],prop_bounds[3],prop_bounds[5]};
			double diagonal = arma::norm(bbox_max - bbox_min);



			if (this -> pos_slider -> value() != 0){
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value() - 1,pos_before.colptr(0));
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value(),pos_now.colptr(0));
			}
			else{
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value(),pos_before.colptr(0));
				wrapped_trajectory_data[trajectory_name] -> get_points() -> GetPoint (this -> pos_slider -> value() + 1,pos_now.colptr(0));

			}

			arma::vec cam_pos = prop_pos_arma + 4 * diagonal * ( arma::normalise(prop_pos_arma) + arma::normalise(pos_before - pos_now));
			camera -> SetPosition(cam_pos.colptr(0));

		}

		this -> parent -> get_renderer() -> Modified();

	}






	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
}




void MoveAlongTrajectoryWindow::prop_removed_slot(){

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_trajectory_data = this -> parent -> get_wrapped_trajectory_data();

	for (unsigned int i = 0; i < this -> spacecraft_combo_box -> count(); ++i){
		QString spacecraft_name = this -> spacecraft_combo_box -> itemText(i);
		if (wrapped_spacecraft_data.find(spacecraft_name.toStdString()) == wrapped_spacecraft_data.end()){
			this -> spacecraft_combo_box -> removeItem(i);
			break;
			// This break is mandatory as the count of spacecraft_combo_box has changed
		}
		
	}
	for (unsigned int i = 0; i < this -> trajectory_combo_box -> count(); ++i){
		QString trajectory_name = this -> trajectory_combo_box -> itemText(i);

		if (wrapped_trajectory_data.find(trajectory_name.toStdString()) == wrapped_trajectory_data.end()){
			this -> trajectory_combo_box -> removeItem(i);
			break;
			// This break is mandatory as the count of trajectory_combo_box has changed

		}
	}

	// If there are no trajectories or spacecraft to use, the slider widget and
	// its spinbox are disabled
	if (this -> trajectory_combo_box -> count() == 0 || this -> spacecraft_combo_box -> count() == 0){
		this -> position_group -> setEnabled(false);
	}



}


void MoveAlongTrajectoryWindow::prop_added_slot(){
	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_trajectory_data = this -> parent -> get_wrapped_trajectory_data();

	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		if (this -> spacecraft_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> spacecraft_combo_box -> insertItem(this -> spacecraft_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}

	for (auto it = wrapped_trajectory_data.begin(); it != wrapped_trajectory_data.end(); ++it){
		if (this -> trajectory_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> trajectory_combo_box -> insertItem(this -> trajectory_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}


	// If there are trajectories and spacecraft to use, the slider widget and
	// its spinbox are enabled
	if (this -> trajectory_combo_box -> count() > 0 && this -> spacecraft_combo_box -> count() > 0){
		this -> position_group -> setEnabled(true);
	}

}

