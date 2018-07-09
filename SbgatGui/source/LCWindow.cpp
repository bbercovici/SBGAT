/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include "LCWindow.hpp"
#include "LCVisualizer.hpp"
#include "ShapePropertiesWidget.hpp"
#include <SBGATMassProperties.hpp>
#include <QMessageBox>
#include <SBGATTrajectory.hpp>
#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>


using namespace SBGAT_GUI;

LCWindow::LCWindow(Mainwindow * parent) : ObsWindow(parent){


	this -> setWindowTitle("Generate Light Curve");

	this -> collect_observations_button = new QPushButton("Collect observations",this);

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	if (wrapped_shape_data.size() == 0){
		QMessageBox::warning(this, "Compute Light Curve", "You must first load a shape model in order to compute light curves!");
		this -> collect_observations_button -> setDisabled(1);
	}

	else{
		this -> collect_observations_button -> setEnabled(1);
	}


	this -> obs_specific_group -> setTitle("Light curve");

	QVBoxLayout * lc_window_layout = new QVBoxLayout(this -> obs_specific_group);

	QGroupBox * observer_location_group = new QGroupBox(tr("Observer location"));
	QGroupBox * sun_location_group = new QGroupBox(tr("Sun location"));
	QGroupBox * phasing_group = new QGroupBox(tr("Phasing"));


	QGridLayout * observer_location_group_layout = new QGridLayout(observer_location_group);
	QGridLayout * sun_location_group_layout = new QGridLayout(sun_location_group);
	QGridLayout * phasing_group_layout = new QGridLayout(phasing_group);


	QLabel * observer_az_label = new QLabel("Observer azimuth (deg)",this);
	QLabel * observer_el_label = new QLabel("Observer elevation (deg)",this);

	QLabel * sun_az_label = new QLabel("Sun azimuth (deg)",this);
	QLabel * sun_el_label = new QLabel("Sun elevation (deg)",this);

	this -> phase_angle_label = new QLabel("Phase angle (deg)" , this);
	this -> phase_angle_qldt = new QLineEdit(this);
	this -> phase_angle_qldt -> setReadOnly(true);
	this -> phase_angle_qldt -> setDisabled(true);

	this -> observer_az_sbox = new QDoubleSpinBox(this);
	this -> observer_el_sbox = new QDoubleSpinBox(this);

	this -> sun_az_sbox = new QDoubleSpinBox(this);
	this -> sun_el_sbox = new QDoubleSpinBox(this);

	observer_location_group_layout -> addWidget(observer_az_label,0,0,1,1);
	observer_location_group_layout -> addWidget(this -> observer_az_sbox,0,1,1,1);

	observer_location_group_layout -> addWidget(observer_el_label,1,0,1,1);
	observer_location_group_layout -> addWidget(this -> observer_el_sbox,1,1,1,1);

	sun_location_group_layout -> addWidget(sun_az_label,0,0,1,1);
	sun_location_group_layout -> addWidget(this -> sun_az_sbox,0,1,1,1);

	sun_location_group_layout -> addWidget(sun_el_label,1,0,1,1);
	sun_location_group_layout -> addWidget(this -> sun_el_sbox,1,1,1,1);

	phasing_group_layout -> addWidget(this -> phase_angle_label,0,0,1,1);
	phasing_group_layout -> addWidget(this -> phase_angle_qldt,0,1,1,1);

	lc_window_layout -> addWidget(observer_location_group);
	lc_window_layout -> addWidget(sun_location_group);
	lc_window_layout -> addWidget(phasing_group);

	lc_window_layout -> addWidget(this -> collect_observations_button);


	connect(this -> collect_observations_button, SIGNAL(clicked()), this, SLOT(collect_observations()));
	connect(this -> open_visualizer_button, SIGNAL(clicked()), this, SLOT(open_visualizer()));
	connect(this -> save_observations_button,SIGNAL(clicked()),this,SLOT(save_observations()));
	connect(this -> observer_az_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> observer_el_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> sun_az_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> sun_el_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));

	this -> init();

}

void LCWindow::init(){

	this -> observer_az_sbox -> setDecimals(6);
	this -> observer_el_sbox -> setDecimals(6);

	this -> sun_az_sbox -> setDecimals(6);
	this -> sun_el_sbox -> setDecimals(6);

	this -> observer_az_sbox -> setRange(-180,180);
	this -> observer_el_sbox -> setRange(-180,180);

	this -> sun_az_sbox -> setRange(-180,180);
	this -> sun_el_sbox -> setRange(-180,180);

	this -> observer_az_sbox -> setValue(0);
	this -> observer_el_sbox -> setValue(0);

	this -> sun_az_sbox -> setValue(0);
	this -> sun_el_sbox -> setValue(0);

	this ->  lc = vtkSmartPointer<SBGATObsLightcurve>::New();
	this -> phase_angle_qldt -> setText(QString::number(0));
	this -> phase_angle_qldt -> repaint();


}



void LCWindow::collect_observations(){


	this -> measurements.clear();


	double imaging_period = this -> imaging_period_sbox -> value() * 3600; 
	std::vector<double> imaging_times;
	for (int i  = 0; i < this -> N_images_sbox -> value(); ++i){
		double t = i * imaging_period;
		imaging_times.push_back(t);
	}




	// Querying the selected primary small body
	std::string primary_name = this -> primary_prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();
	this -> lc -> SetInputData(0,shape_data[primary_name] -> get_polydata());
	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	mass_properties -> SetInputData(shape_data[primary_name] -> get_polydata());
	mass_properties -> Update();
	double mu = arma::datum::G * this -> primary_shape_properties_widget -> get_density() * mass_properties -> GetVolume();
	
	// Querying the selected secondary small body, if any
	std::string secondary_name = this -> secondary_prop_combo_box -> currentText().toStdString();
	arma::vec elements;
	SBGATTrajectory trajectory;
	std::vector<arma::vec> secondary_positions,secondary_velocities;

	if (secondary_name != "None"){
		shape_data = this -> parent -> get_wrapped_shape_data();
		this -> lc -> SetInputData(1,shape_data[secondary_name] -> get_polydata());
		elements = this -> secondary_shape_properties_widget -> get_orbital_elements();
		trajectory.GenerateKeplerianTrajectory(secondary_positions,secondary_velocities,imaging_times,elements,mu);
	}


	this -> lc -> SetScaleMeters();
	this -> lc -> Update();


	double d2r = arma::datum::pi /180;
	arma::vec observer_dir = {1,0,0};
	observer_dir = (RBK::M2(this -> observer_el_sbox -> value() * d2r) 
		* RBK::M3(this -> observer_az_sbox -> value() * d2r)).t() * observer_dir;
	arma::vec sun_dir = {1,0,0};
	sun_dir = (RBK::M2(this -> sun_el_sbox -> value() * d2r) 
		* RBK::M3(this -> sun_el_sbox -> value() * d2r)).t() * sun_dir;


	for (unsigned int t = 0; t < imaging_times.size(); ++t){

		std::vector<double> period_vec;
		std::vector<arma::vec> positions_vec, velocities_vec, spin_vec;

		// Primary
		period_vec.push_back(this -> primary_shape_properties_widget -> get_period());
		positions_vec.push_back(arma::zeros<arma::vec>(3));
		velocities_vec.push_back(arma::zeros<arma::vec>(3));
		spin_vec.push_back(this -> primary_shape_properties_widget -> get_spin());

		// Secondary, if any
		if (secondary_name != "None"){
			period_vec.push_back(this -> secondary_shape_properties_widget -> get_period());
			positions_vec.push_back(secondary_positions[t]);
			velocities_vec.push_back(secondary_velocities[t]);
			spin_vec.push_back(this -> secondary_shape_properties_widget -> get_spin());
		}

		this -> lc -> CollectMeasurementsSimpleSpin(
			this -> measurements,
			this -> N_samples_sbox -> value(),
			imaging_times[t],
			period_vec,
			sun_dir,
			observer_dir,
			positions_vec,
			spin_vec,
			this -> penalize_incidence_box -> isChecked());
	}


	this -> open_visualizer_button -> setEnabled(1);
	this -> save_observations_button -> setEnabled(1);

	this -> open_visualizer_button -> repaint();
	this -> save_observations_button -> repaint();

}

void LCWindow::open_visualizer(){

	LCVisualizer lc_visualizer(this,this -> measurements);
	lc_visualizer.exec();

}

void LCWindow::save_observations(){

	QString path = QFileDialog::getSaveFileName(this, tr("Save File"),
		"",
		tr("Text file (*.txt)"));
	if (path.size() != 0){
		this -> lc -> SaveLightCurveData(this -> measurements, path.toStdString());
	}

}

void LCWindow::update_phase_angle(){

	double d2r = arma::datum::pi /180;
	double r2d = 1./d2r;

	
	arma::vec observer_dir = {1,0,0};
	observer_dir = (RBK::M2(this -> observer_el_sbox -> value() * d2r) 
		* RBK::M3(this -> observer_az_sbox -> value() * d2r)).t() * observer_dir;

	arma::vec sun_dir = {1,0,0};
	sun_dir = (RBK::M2(this -> sun_el_sbox -> value() * d2r) 
		* RBK::M3(this -> sun_az_sbox -> value() * d2r)).t() * sun_dir;

	double phase_angle = std::acos(arma::dot(observer_dir,sun_dir)) * r2d;

	this -> phase_angle_qldt -> setText(QString::number(phase_angle));
	this -> phase_angle_qldt -> repaint();

}	


