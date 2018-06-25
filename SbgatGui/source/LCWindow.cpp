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
#include <QMessageBox>


using namespace SBGAT_GUI;

LCWindow::LCWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Generate Light Curve");

	this -> save_observations_button = new QPushButton("Save observations",this);
	this -> collect_observations_button = new QPushButton("Collect observations",this);
	
	this -> open_visualizer_button = new QPushButton("Visualize observations",this);


	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	if (wrapped_shape_data.size() == 0){
		QMessageBox::warning(this, "Compute Light Curve", "You must first load a shape model in order to compute light curves!");
		collect_observations_button -> setDisabled(1);
	}

	else{
		collect_observations_button -> setEnabled(1);
	}


	


	QVBoxLayout * lc_window_layout = new QVBoxLayout(this);

	QGroupBox * target_group = new QGroupBox(tr("Shape"));
	QGroupBox * target_properties_group = new QGroupBox(tr("Shape properties"));
	QGroupBox * lc_settings_group = new QGroupBox(tr("Light curve settings"));

	QGridLayout * target_group_layout = new QGridLayout(target_group);
	QGridLayout * lc_settings_group_layout = new QGridLayout(lc_settings_group);
	QGridLayout * target_properties_group_layout = new QGridLayout(target_properties_group);

	QLabel * shape_label = new QLabel("Targeted shape model",this);
	QLabel * N_samples_label = new QLabel("Max samples per facet",this);

	QLabel * spin_raan_label = new QLabel("Spin RAAN (deg)",this);
	QLabel * spin_inc_label = new QLabel("Spin inclination (deg)",this);

	QLabel * observer_az_label = new QLabel("Observer azimuth (deg)",this);
	QLabel * observer_el_label = new QLabel("Observer elevation (deg)",this);

	QLabel * sun_az_label = new QLabel("Sun azimuth (deg)",this);
	QLabel * sun_el_label = new QLabel("Sun elevation (deg)",this);

	QLabel * period_label = new QLabel("Rotation period (hours)",this);
	QLabel * imaging_period_label = new QLabel("Imaging period (hours)",this);
	QLabel * N_images_label = new QLabel("Images to collect",this);


	this -> phase_angle_label = new QLabel("Phase angle (deg)" , this);
	this -> phase_angle_qldt = new QLineEdit(this);
	this -> phase_angle_qldt -> setReadOnly(true);
	this -> phase_angle_qldt -> setDisabled(true);

	this -> prop_combo_box = new QComboBox (this);
	this -> spin_raan_sbox = new QDoubleSpinBox(this);
	this -> spin_inc_sbox = new QDoubleSpinBox(this);

	this -> observer_az_sbox = new QDoubleSpinBox(this);
	this -> observer_el_sbox = new QDoubleSpinBox(this);


	this -> sun_az_sbox = new QDoubleSpinBox(this);
	this -> sun_el_sbox = new QDoubleSpinBox(this);

	this -> rotation_period_sbox = new QDoubleSpinBox(this);
	this -> imaging_period_sbox = new QDoubleSpinBox(this);

	this -> N_samples_sbox = new QSpinBox (this);
	this -> N_images_sbox = new QSpinBox (this);

	this -> penalize_incidence_box = new QCheckBox("Penalize incidence",this);


	target_group_layout -> addWidget(shape_label,0,0,1,1);
	target_group_layout -> addWidget(this -> prop_combo_box,0,1,1,1);

	target_properties_group_layout -> addWidget(spin_raan_label,0,0,1,1);
	target_properties_group_layout -> addWidget(this -> spin_raan_sbox,0,1,1,1);

	target_properties_group_layout -> addWidget(spin_inc_label,1,0,1,1);
	target_properties_group_layout -> addWidget(this -> spin_inc_sbox,1,1,1,1);

	target_properties_group_layout -> addWidget(period_label,2,0,1,1);
	target_properties_group_layout -> addWidget(this -> rotation_period_sbox,2,1,1,1);


	lc_settings_group_layout -> addWidget(N_samples_label,0,0,1,1);
	lc_settings_group_layout -> addWidget(this -> N_samples_sbox,0,1,1,1);

	lc_settings_group_layout -> addWidget(imaging_period_label,1,0,1,1);
	lc_settings_group_layout -> addWidget(this -> imaging_period_sbox,1,1,1,1);

	lc_settings_group_layout -> addWidget(N_images_label,2,0,1,1);
	lc_settings_group_layout -> addWidget(this -> N_images_sbox,2,1,1,1);

	lc_settings_group_layout -> addWidget(observer_az_label,3,0,1,1);
	lc_settings_group_layout -> addWidget(this -> observer_az_sbox,3,1,1,1);

	lc_settings_group_layout -> addWidget(observer_el_label,4,0,1,1);
	lc_settings_group_layout -> addWidget(this -> observer_el_sbox,4,1,1,1);

	lc_settings_group_layout -> addWidget(sun_az_label,5,0,1,1);
	lc_settings_group_layout -> addWidget(this -> sun_az_sbox,5,1,1,1);

	lc_settings_group_layout -> addWidget(sun_el_label,6,0,1,1);
	lc_settings_group_layout -> addWidget(this -> sun_el_sbox,6,1,1,1);

	lc_settings_group_layout -> addWidget(this -> phase_angle_label,7,0,1,1);
	lc_settings_group_layout -> addWidget(this -> phase_angle_qldt,7,1,1,1);
	lc_settings_group_layout -> addWidget(this -> penalize_incidence_box,8,0,1,2);


	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	lc_window_layout -> addWidget(target_group);
	lc_window_layout -> addWidget(target_properties_group);
	lc_window_layout -> addWidget(lc_settings_group);
	lc_window_layout -> addWidget(collect_observations_button);
	lc_window_layout -> addWidget(open_visualizer_button);
	lc_window_layout -> addWidget(save_observations_button);


	open_visualizer_button -> setDisabled(1);
	save_observations_button -> setDisabled(1);

	lc_window_layout -> addWidget(button_box);

	connect(collect_observations_button, SIGNAL(clicked()), this, SLOT(collect_observations()));
	connect(open_visualizer_button, SIGNAL(clicked()), this, SLOT(open_visualizer()));
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> save_observations_button,SIGNAL(clicked()),this,SLOT(save_observations()));
	connect(this -> observer_az_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> observer_el_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> sun_az_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));
	connect(this -> sun_el_sbox,SIGNAL(valueChanged(double)),this,SLOT(update_phase_angle()));


	this -> init();


}

void LCWindow::init(){

	this -> spin_raan_sbox -> setDecimals(6);
	this -> spin_inc_sbox -> setDecimals(6);

	this -> observer_az_sbox -> setDecimals(6);
	this -> observer_el_sbox -> setDecimals(6);

	this -> sun_az_sbox -> setDecimals(6);
	this -> sun_el_sbox -> setDecimals(6);

	this -> rotation_period_sbox -> setDecimals(6);
	this -> imaging_period_sbox -> setDecimals(6);

	this -> spin_raan_sbox -> setRange(-180,180);
	this -> spin_inc_sbox -> setRange(-180,180);

	this -> observer_az_sbox -> setRange(-180,180);
	this -> observer_el_sbox -> setRange(-180,180);


	this -> sun_az_sbox -> setRange(-180,180);
	this -> sun_el_sbox -> setRange(-180,180);

	this -> N_samples_sbox -> setRange(1,1000);
	this -> N_images_sbox -> setRange(1,1000);

	this -> rotation_period_sbox -> setRange(1e-10,1e10);
	this -> imaging_period_sbox -> setRange(1e-10,1e10);

	this -> spin_raan_sbox -> setValue(0);
	this -> spin_inc_sbox -> setValue(0);

	this -> observer_az_sbox -> setValue(0);
	this -> observer_el_sbox -> setValue(0);


	this -> sun_az_sbox -> setValue(0);
	this -> sun_el_sbox -> setValue(0);

	this -> N_samples_sbox -> setValue(30);
	this -> N_images_sbox -> setValue(1);

	this -> rotation_period_sbox -> setValue(1);
	this -> imaging_period_sbox -> setValue(1);

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
	}

	if (wrapped_shape_data.size() == 0 ){
		this -> save_observations_button -> setEnabled(false);
	}

	this ->  lc = vtkSmartPointer<SBGATObsLightcurve>::New();
	this -> phase_angle_qldt -> setText(QString::number(0));
	this -> phase_angle_qldt -> repaint();


}



void LCWindow::collect_observations(){

	double d2r = arma::datum::pi /180;
	arma::vec spin = {0,0,1};
	spin = (RBK::M2(this -> spin_inc_sbox -> value() * d2r) 
		* RBK::M3(this -> spin_raan_sbox -> value() * d2r)).t() * spin;
	
	arma::vec observer_dir = {1,0,0};
	observer_dir = (RBK::M2(this -> observer_el_sbox -> value() * d2r) 
		* RBK::M3(this -> observer_az_sbox -> value() * d2r)).t() * observer_dir;

	arma::vec sun_dir = {1,0,0};
	sun_dir = (RBK::M2(this -> sun_el_sbox -> value() * d2r) 
		* RBK::M3(this -> sun_el_sbox -> value() * d2r)).t() * sun_dir;

	double rotation_period = this -> rotation_period_sbox -> value() * 3600; 
	double imaging_period = this -> imaging_period_sbox -> value() * 3600; 

	int N_images = this -> N_images_sbox -> value(); 
	int N_samples = this -> N_samples_sbox -> value() ;


	// Querying the selected small body
	std::string name = this -> prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();

	this -> lc -> SetInputData(shape_data[name] -> get_polydata());
	this -> lc -> SetScaleMeters();
	this -> lc -> Update();

	this -> measurements.clear();

	for (int i  = 0; i < N_images; ++i){

		double t = i * imaging_period;
		this -> lc -> CollectMeasurementsSimpleSpin(this -> measurements,
			N_samples,
			t,
			rotation_period,
			sun_dir,
			observer_dir,
			spin,
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


