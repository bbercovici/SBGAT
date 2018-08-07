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

	this -> phase_angle_qldt -> setText(QString::number(0));
	this -> phase_angle_qldt -> repaint();


}



void LCWindow::collect_observations(){


	this -> measurements.clear();

	this -> observation_filter = vtkSmartPointer<SBGATObsLightcurve>::New();

	std::vector<double> imaging_times;
	std::vector< std::vector<arma::vec> > positions_vec, velocities_vec, mrps_vec,omegas_vec;

	ObsWindow::get_inputs_from_GUI(imaging_times,
		positions_vec,
		velocities_vec,
		mrps_vec,
		omegas_vec);
	
	this -> observation_filter -> SetScaleMeters();

	double d2r = arma::datum::pi /180;
	arma::vec observer_dir = {1,0,0};
	observer_dir = (RBK::M2(this -> observer_el_sbox -> value() * d2r) 
		* RBK::M3(this -> observer_az_sbox -> value() * d2r)).t() * observer_dir;
	arma::vec sun_dir = {1,0,0};
	sun_dir = (RBK::M2(this -> sun_el_sbox -> value() * d2r) 
		* RBK::M3(this -> sun_el_sbox -> value() * d2r)).t() * sun_dir;

	SBGATObsLightcurve * lc = SBGATObsLightcurve::SafeDownCast(this -> observation_filter);


	for (unsigned int t = 0; t < imaging_times.size(); ++t){

		
		lc -> CollectMeasurements(
			this -> measurements,
			imaging_times[t],
			this -> N_samples_sbox -> value(),
			sun_dir,
			observer_dir,
			positions_vec[t],
			velocities_vec[t], 
			mrps_vec[t],
			omegas_vec[t],
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
	SBGATObsLightcurve * lc = SBGATObsLightcurve::SafeDownCast(this -> observation_filter);
	
	QString path = QFileDialog::getSaveFileName(this, tr("Save File"),
		"",
		tr("Text file (*.txt)"));
	if (path.size() != 0){
		lc -> SaveLightCurveData(this -> measurements, path.toStdString());
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


