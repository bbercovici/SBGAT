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

#include "RadarWindow.hpp"
#include <SBGATObsRadar.hpp>

using namespace SBGAT_GUI;

RadarWindow::RadarWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Generate Doppler Radar Observations");


	QVBoxLayout * radar_window_layout = new QVBoxLayout(this);

	QGroupBox * target_group = new QGroupBox(tr("Shape"));
	QGroupBox * target_properties_group = new QGroupBox(tr("Shape properties"));
	QGroupBox * radar_settings_group = new QGroupBox(tr("Radar settings"));
	QGroupBox * output_group = new QGroupBox(tr("Radar output"));
	
	QGridLayout * target_group_layout = new QGridLayout(target_group);
	QGridLayout * radar_settings_group_layout = new QGridLayout(radar_settings_group);
	QGridLayout * target_properties_group_layout = new QGridLayout(target_properties_group);
	QGridLayout * output_group_layout = new QGridLayout(output_group);

	QLabel * shape_label = new QLabel("Targeted shape model",this);
	QLabel * range_label = new QLabel("Range resolution (m)",this);
	QLabel * range_rate_label = new QLabel("Range-rate resolution (m/s)",this);
	QLabel * N_samples_label = new QLabel("Max samples per facet",this);


	QLabel * spin_az_label = new QLabel("Spin azimuth (deg)",this);
	QLabel * spin_el_label = new QLabel("Spin elevation (deg)",this);
	QLabel * period_label = new QLabel("Rotation period (hours)",this);
	QLabel * imaging_period_label = new QLabel("Imaging period (hours)",this);
	QLabel * N_images_label = new QLabel("Images to collect",this);



	this -> prop_combo_box = new QComboBox (this);
	this -> r_bin_sbox = new QDoubleSpinBox(this);
	this -> rr_bin_sbox = new QDoubleSpinBox(this);
	this -> spin_az_sbox = new QDoubleSpinBox(this);
	this -> spin_el_sbox = new QDoubleSpinBox(this);
	this -> rotation_period_sbox = new QDoubleSpinBox(this);
	this -> imaging_period_sbox = new QDoubleSpinBox(this);

	this -> N_samples_sbox = new QSpinBox (this);
	this -> N_images_sbox = new QSpinBox (this);

	this ->  open_output_file_dialog_button = new QPushButton("Select output folder",this);
	this -> collect_observations_button = new QPushButton("Collect observations",this);

	target_group_layout -> addWidget(shape_label,0,0,1,1);
	target_group_layout -> addWidget(this -> prop_combo_box,0,1,1,1);

	target_properties_group_layout -> addWidget(spin_az_label,0,0,1,1);
	target_properties_group_layout -> addWidget(this -> spin_az_sbox,0,1,1,1);

	target_properties_group_layout -> addWidget(spin_el_label,1,0,1,1);
	target_properties_group_layout -> addWidget(this -> spin_el_sbox,1,1,1,1);

	target_properties_group_layout -> addWidget(period_label,2,0,1,1);
	target_properties_group_layout -> addWidget(this -> rotation_period_sbox,2,1,1,1);

	radar_settings_group_layout -> addWidget(range_label,0,0,1,1);
	radar_settings_group_layout -> addWidget(this -> r_bin_sbox,0,1,1,1);

	radar_settings_group_layout -> addWidget(range_rate_label,1,0,1,1);
	radar_settings_group_layout -> addWidget(this -> rr_bin_sbox,1,1,1,1);

	radar_settings_group_layout -> addWidget(N_samples_label,2,0,1,1);
	radar_settings_group_layout -> addWidget(this -> N_samples_sbox,2,1,1,1);

	radar_settings_group_layout -> addWidget(imaging_period_label,3,0,1,1);
	radar_settings_group_layout -> addWidget(this -> imaging_period_sbox,3,1,1,1);


	radar_settings_group_layout -> addWidget(N_images_label,4,0,1,1);
	radar_settings_group_layout -> addWidget(this -> N_images_sbox,4,1,1,1);
	
	output_group_layout -> addWidget(this -> open_output_file_dialog_button, 0, 0, 1, 2);


	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	radar_window_layout -> addWidget(target_group);
	radar_window_layout -> addWidget(target_properties_group);
	radar_window_layout -> addWidget(radar_settings_group);
	radar_window_layout -> addWidget(output_group);
	radar_window_layout -> addWidget(collect_observations_button);
	radar_window_layout -> addWidget(button_box);


	connect(collect_observations_button, SIGNAL(clicked()), this, SLOT(collect_observations()));

	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> open_output_file_dialog_button,SIGNAL(clicked()),this,
		SLOT(open_output_file_dialog()));

	this -> init();


}

void RadarWindow::init(){


	this -> r_bin_sbox -> setDecimals(6);
	this -> rr_bin_sbox -> setDecimals(6);
	this -> spin_az_sbox -> setDecimals(6);
	this -> spin_el_sbox -> setDecimals(6);
	this -> rotation_period_sbox -> setDecimals(6);
	this -> imaging_period_sbox -> setDecimals(6);

	this -> r_bin_sbox -> setRange(1e-10,1e10);
	this -> rr_bin_sbox -> setRange(1e-10,1e10);
	this -> spin_az_sbox -> setRange(-180,180);
	this -> spin_el_sbox -> setRange(-180,180);
	this -> N_samples_sbox -> setRange(1,1000);
	this -> N_images_sbox -> setRange(1,1000);

	this -> rotation_period_sbox -> setRange(1e-10,1e10);
	this -> imaging_period_sbox -> setRange(1e-10,1e10);



	this -> r_bin_sbox -> setValue(2);
	this -> rr_bin_sbox -> setValue(1e-3);
	this -> spin_az_sbox -> setValue(0);
	this -> spin_el_sbox -> setValue(0);
	this -> N_samples_sbox -> setValue(100);
	this -> N_images_sbox -> setValue(10);

	this -> rotation_period_sbox -> setValue(1);
	this -> imaging_period_sbox -> setValue(1);



	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
	}

	
	this -> collect_observations_button -> setEnabled(false);

	if (wrapped_shape_data.size() == 0 ){
		this -> open_output_file_dialog_button -> setEnabled(false);
	}


}

void RadarWindow::open_output_file_dialog(){

	std::string default_name;

	std::string name = this -> prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();
	
	if ( shape_data.find(name)!= shape_data.end()){
		default_name = name;
	}

	QString path = QFileDialog::getExistingDirectory(this, tr("Select radar output folder"),
		QString::fromStdString(default_name));


	this -> output_path = path.toStdString();

	if (this -> output_path.size() > 0){
		this -> collect_observations_button -> setEnabled(true);
	}
}

void RadarWindow::collect_observations(){

	double d2r = arma::datum::pi /180;
	arma::vec spin = {0,0,1};
	spin = (RBK::M2(this -> spin_el_sbox -> value() * d2r) 
		* RBK::M3(this -> spin_az_sbox -> value() * d2r)).t() * spin;
	
	arma::vec dir = {-1,0,0};

	double rotation_period = this -> rotation_period_sbox -> value() * 3600; 
	double imaging_period = this -> imaging_period_sbox -> value() * 3600; 

	
	int N_images = this -> N_images_sbox -> value(); 
	int N_samples = this -> N_samples_sbox -> value() ;

	double r_bin = this -> r_bin_sbox -> value();
	double rr_bin = this -> rr_bin_sbox -> value();

	// Querying the selected small body
	std::string name = this -> prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();

	vtkSmartPointer<SBGATObsRadar> radar = vtkSmartPointer<SBGATObsRadar>::New();
	radar -> SetInputData(shape_data[name] -> get_polydata());
	radar -> SetScaleMeters();
	radar -> Update();


	for (int i  = 0; i < N_images; ++i){

		double t = i * imaging_period;

		std::vector<std::array<double, 2> > measurements;

		radar -> CollectMeasurementsSimpleSpin(measurements,N_samples,t,rotation_period,dir,spin);
		radar -> BinObs(measurements,r_bin,rr_bin);
		radar -> SaveImage( this -> output_path + "/" + std::to_string(i) +  ".png");

	}

}



