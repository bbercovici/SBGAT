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
#include "RadarVisualizer.hpp"
#include "ShapePropertiesWidget.hpp"

#include <QMessageBox>

#include <SBGATMassProperties.hpp>



using namespace SBGAT_GUI;

RadarWindow::RadarWindow(Mainwindow * parent) : ObsWindow(parent) {


	this -> parent = parent;
	this -> setWindowTitle("Generate Doppler Radar Observations");

	this -> collect_observations_button = new QPushButton("Collect observations",this);
	this -> bin_observations_button = new QPushButton("Bin observations",this);

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	collect_observations_button -> setEnabled(1);
	


	this -> obs_specific_group -> setTitle("Radar");

	QVBoxLayout * radar_window_layout = new QVBoxLayout(this -> obs_specific_group);

	QGroupBox * radar_settings_group = new QGroupBox(tr("Radar location"));
	QGroupBox * binning_settings_group = new QGroupBox(tr("Binning settings"));
	
	QGridLayout * radar_settings_group_layout = new QGridLayout(radar_settings_group);
	QGridLayout * binning_settings_group_layout = new QGridLayout(binning_settings_group);


	QLabel * range_label = new QLabel("Range resolution (m)",this);
	QLabel * range_rate_label = new QLabel("Range-rate resolution (m/s)",this);

	QLabel * radar_az_label = new QLabel("Radar azimuth (deg)",this);
	QLabel * radar_el_label = new QLabel("Radar elevation (deg)",this);

	
	this -> r_bin_sbox = new QDoubleSpinBox(this);
	this -> rr_bin_sbox = new QDoubleSpinBox(this);
	
	this -> radar_az_sbox = new QDoubleSpinBox(this);
	this -> radar_el_sbox = new QDoubleSpinBox(this);

	
	radar_settings_group_layout -> addWidget(radar_az_label,3,0,1,1);
	radar_settings_group_layout -> addWidget(this -> radar_az_sbox,3,1,1,1);

	radar_settings_group_layout -> addWidget(radar_el_label,4,0,1,1);
	radar_settings_group_layout -> addWidget(this -> radar_el_sbox,4,1,1,1);


	binning_settings_group_layout -> addWidget(range_label,0,0,1,1);
	binning_settings_group_layout -> addWidget(this -> r_bin_sbox,0,1,1,1);

	binning_settings_group_layout -> addWidget(range_rate_label,1,0,1,1);
	binning_settings_group_layout -> addWidget(this -> rr_bin_sbox,1,1,1,1);
	

	radar_window_layout -> addWidget(radar_settings_group);
	radar_window_layout -> addWidget(collect_observations_button);
	radar_window_layout -> addWidget(binning_settings_group);
	radar_window_layout -> addWidget(bin_observations_button);

	bin_observations_button -> setDisabled(1);

	connect(collect_observations_button, SIGNAL(clicked()), this, SLOT(collect_observations()));
	connect(bin_observations_button, SIGNAL(clicked()), this, SLOT(bin_observations()));
	

	connect(this -> open_visualizer_button, SIGNAL(clicked()), this, SLOT(open_visualizer()));
	connect(this -> save_observations_button,SIGNAL(clicked()),this,SLOT(save_observations()));


	this -> init();


}

void RadarWindow::init(){

	this -> r_bin_sbox -> setDecimals(6);
	this -> rr_bin_sbox -> setDecimals(6);
	
	this -> radar_az_sbox -> setDecimals(6);
	this -> radar_el_sbox -> setDecimals(6);
	
	this -> r_bin_sbox -> setRange(1e-10,1e10);
	this -> rr_bin_sbox -> setRange(1e-10,1e10);
	
	this -> radar_az_sbox -> setRange(-180,180);
	this -> radar_el_sbox -> setRange(-180,180);

	this -> r_bin_sbox -> setValue(2);
	this -> rr_bin_sbox -> setValue(1e-3);
	

	this -> radar_az_sbox -> setValue(0);
	this -> radar_el_sbox -> setValue(0);
	
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	if (wrapped_shape_data.size() == 0 ){
		this -> save_observations_button -> setEnabled(false);
	}


}


void RadarWindow::open_visualizer(){

	SBGATObsRadar * radar = SBGATObsRadar::SafeDownCast(this -> observation_filter);

	RadarVisualizer radar_visualizer(this,radar -> GetImages());
	radar_visualizer.exec();

}


void RadarWindow::collect_observations(){
	
	this -> measurement_sequence.clear();

	this -> observation_filter = vtkSmartPointer<SBGATObsRadar>::New();
	
	std::vector<double> imaging_times;
	std::vector< std::vector<arma::vec> > positions_vec, velocities_vec, mrps_vec,omegas_vec;

	this -> get_inputs_from_GUI(imaging_times,
		positions_vec,
		velocities_vec,
		mrps_vec,
		omegas_vec);

	this -> observation_filter -> SetScaleMeters();
	this -> observation_filter -> Update();


	double d2r = arma::datum::pi /180;
	arma::vec radar_dir = {1,0,0};
	radar_dir = (RBK::M2(this -> radar_el_sbox -> value() * d2r) 
		* RBK::M3(this -> radar_az_sbox -> value() * d2r)).t() * radar_dir;


	SBGATObsRadar * radar = SBGATObsRadar::SafeDownCast(this -> observation_filter);

	for (unsigned int t = 0; t < imaging_times.size(); ++t){

		radar -> CollectMeasurements(this -> measurement_sequence,
			imaging_times[t],
			this -> N_samples_sbox -> value(),
			radar_dir,
			positions_vec[t],
			velocities_vec[t], 
			mrps_vec[t],
			omegas_vec[t],
			this -> penalize_incidence_box -> isChecked());
	}

	

	this -> bin_observations_button -> setEnabled(1);
	this -> bin_observations_button -> repaint();

}


void RadarWindow::bin_observations(){


	double r_bin = this -> r_bin_sbox -> value();
	double rr_bin = this -> rr_bin_sbox -> value();
	SBGATObsRadar * radar = SBGATObsRadar::SafeDownCast(this -> observation_filter);

	if (r_bin == 0 || rr_bin == 0){


		radar -> ClearImages();

		QMessageBox::warning(this, "Generate Doppler Radar Observations", "Invalid bin size");
		this -> open_visualizer_button -> setEnabled(0);
		this -> save_observations_button -> setEnabled(0);

		this -> open_visualizer_button -> repaint();
		this -> save_observations_button -> repaint();

	}
	else{
		try{
			radar -> BinObservations(this -> measurement_sequence,r_bin,rr_bin);
			this -> open_visualizer_button -> setEnabled(1);
			this -> save_observations_button -> setEnabled(1);

			this -> open_visualizer_button -> repaint();
			this -> save_observations_button -> repaint();
		}

		catch(std::runtime_error & e){

			QMessageBox::warning(this, "Generate Doppler Radar Observations", e.what());

			this -> open_visualizer_button -> setEnabled(0);
			this -> save_observations_button -> setEnabled(0);

			this -> open_visualizer_button -> repaint();
			this -> save_observations_button -> repaint();


		}
		
		
	}


}

void RadarWindow::save_observations(){
	SBGATObsRadar * radar = SBGATObsRadar::SafeDownCast(this -> observation_filter);

	QString path = QFileDialog::getExistingDirectory(this, tr("Select output folder"));
	if (path.size() != 0){
		radar -> SaveImages( path.toStdString() + "/");
	}

}



