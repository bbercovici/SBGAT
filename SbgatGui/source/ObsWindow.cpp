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

#include <QScrollArea>

#include "ShapePropertiesWidget.hpp"
#include "ObsWindow.hpp"

#include <SBGATMassProperties.hpp>
#include <SBGATTrajectory.hpp>
#include <OrbitConversions.hpp>


using namespace SBGAT_GUI;

ObsWindow::ObsWindow(Mainwindow * parent) {

	this -> parent = parent;


	QGroupBox * target_group = new QGroupBox(tr("Shapes"));

	QScrollArea * scroll_area = new QScrollArea(this);
	scroll_area -> setFrameShape(QFrame::NoFrame);
	scroll_area -> setWidgetResizable( true );
	scroll_area -> setGeometry( 10, 10, 600, 800 );

	QWidget * enclosing_widget = new QWidget();

	QVBoxLayout * obs_window_layout = new QVBoxLayout(enclosing_widget);
	scroll_area -> setWidget(enclosing_widget);
	

	this -> obs_specific_group = new QGroupBox(tr("-"));
	QGroupBox * settings_group = new QGroupBox(tr("Settings"));

	QGridLayout * target_group_layout = new QGridLayout(target_group);
	QGridLayout * settings_group_layout = new QGridLayout(settings_group);
	
	QLabel * primary_shape_label = new QLabel("Primary shape model",this);
	QLabel * secondary_shape_label = new QLabel("Secondary shape model",this);

	QLabel * N_samples_label = new QLabel("Minimum number of samples per facet",this);


	QLabel * imaging_period_label = new QLabel("Imaging period (hours)",this);
	QLabel * N_images_label = new QLabel("Images to collect",this);

	this -> open_visualizer_button = new QPushButton("Visualize observations",this);
	this -> save_observations_button = new QPushButton("Save observations",this);
	
	this -> primary_prop_combo_box = new QComboBox (this);
	this -> secondary_prop_combo_box = new QComboBox (this);

	this -> imaging_period_sbox = new QDoubleSpinBox(this);


	this -> N_samples_sbox = new QSpinBox (this);
	this -> N_images_sbox = new QSpinBox (this);

	this -> penalize_incidence_box = new QCheckBox("Penalize incidence",this);

	target_group_layout -> addWidget(primary_shape_label,0,0,1,1);
	target_group_layout -> addWidget(this -> primary_prop_combo_box,0,1,1,1);

	target_group_layout -> addWidget(secondary_shape_label,1,0,1,1);
	target_group_layout -> addWidget(this -> secondary_prop_combo_box,1,1,1,1);

	settings_group_layout -> addWidget(N_samples_label,0,0,1,1);
	settings_group_layout -> addWidget(this -> N_samples_sbox,0,1,1,1);

	settings_group_layout -> addWidget(imaging_period_label,1,0,1,1);
	settings_group_layout -> addWidget(this -> imaging_period_sbox,1,1,1,1);

	settings_group_layout -> addWidget(N_images_label,2,0,1,1);
	settings_group_layout -> addWidget(this -> N_images_sbox,2,1,1,1);

	settings_group_layout -> addWidget(this -> penalize_incidence_box,3,0,1,2);


	this -> primary_shape_properties_widget = new ShapePropertiesWidget(this ,true,"Primary shape properties");
	
	this -> secondary_shape_properties_widget = new ShapePropertiesWidget(this ,false,"Secondary shape properties");

	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	obs_window_layout -> addWidget(target_group);
	obs_window_layout -> addWidget(primary_shape_properties_widget);
	obs_window_layout -> addWidget(secondary_shape_properties_widget);
	obs_window_layout -> addWidget(settings_group);
	obs_window_layout -> addWidget(this -> obs_specific_group);
	obs_window_layout -> addWidget(open_visualizer_button);
	obs_window_layout -> addWidget(save_observations_button);
	obs_window_layout -> addWidget(button_box);


	
	this -> init();
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(secondary_prop_combo_box,SIGNAL(currentIndexChanged(int)),this,SLOT(changed_secondary_box(int)));


}

void ObsWindow::init(){

	this -> imaging_period_sbox -> setDecimals(6);


	this -> N_samples_sbox -> setRange(1,1000);
	this -> N_images_sbox -> setRange(1,1000);

	this -> imaging_period_sbox -> setRange(1e-10,1e10);


	this -> N_samples_sbox -> setValue(1);
	this -> N_images_sbox -> setValue(1);

	this -> imaging_period_sbox -> setValue(1);

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	this -> secondary_prop_combo_box -> insertItem(0,"None");

	if (wrapped_shape_data.size() == 0 ){
		this -> save_observations_button -> setEnabled(false);
	}
	else{
		for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
			this -> primary_prop_combo_box -> insertItem(this -> primary_prop_combo_box -> count(),QString::fromStdString(it -> first));
			this -> secondary_prop_combo_box -> insertItem(this -> secondary_prop_combo_box -> count(),QString::fromStdString(it -> first));
		}
	}

	this -> penalize_incidence_box -> setChecked(true);
	this -> open_visualizer_button -> setDisabled(true);
	this -> save_observations_button -> setDisabled(true);
	this -> secondary_shape_properties_widget -> setEnabled(0);

}


void ObsWindow::changed_secondary_box(int index){
	if (index != 0){
		this -> secondary_shape_properties_widget -> setEnabled(1);
	}
	else{
		this -> secondary_shape_properties_widget -> setEnabled(0);
	}
}


void ObsWindow::get_inputs_from_GUI(std::vector<double> & imaging_times,
	std::vector< std::vector<arma::vec> > & positions_vec, 
	std::vector< std::vector<arma::vec> > & velocities_vec,
	std::vector< std::vector<arma::vec> > & mrps_vec,
	std::vector< std::vector<arma::vec> > & omegas_vec){

	std::string primary_name = this -> primary_prop_combo_box -> currentText().toStdString();
	std::string secondary_name = this -> secondary_prop_combo_box -> currentText().toStdString();
		
	for (int i  = 0; i < this -> N_images_sbox -> value(); ++i){
		double time = i * this -> imaging_period_sbox -> value() * 3600;
		imaging_times.push_back(time);
		positions_vec.push_back(std::vector<arma::vec>());
		velocities_vec.push_back(std::vector<arma::vec>());
		mrps_vec.push_back(std::vector<arma::vec>());
		omegas_vec.push_back(std::vector<arma::vec>());
	}


	// Primary state
	auto shape_data = this -> parent -> get_wrapped_shape_data();

	this -> observation_filter -> AddInputData(0,shape_data[primary_name] -> get_polydata());

	this -> add_state_history(this -> primary_shape_properties_widget,

		imaging_times,positions_vec,velocities_vec,mrps_vec,omegas_vec);


	// If there's a secondary, its state history is  added
	if (secondary_name != "None"){
		this -> observation_filter -> AddInputData(0,shape_data[secondary_name] -> get_polydata());
		this -> add_state_history(this -> secondary_shape_properties_widget,
			imaging_times,positions_vec, velocities_vec,mrps_vec,omegas_vec);
	}




}


void ObsWindow::add_state_history(ShapePropertiesWidget * shape_properties_widget,
	const std::vector<double> & imaging_times,
	std::vector< std::vector<arma::vec> > & positions_vec, 
	std::vector< std::vector<arma::vec> > & velocities_vec,
	std::vector< std::vector<arma::vec> > & mrps_vec,
	std::vector< std::vector<arma::vec> > & omegas_vec) const{


	auto shape_data = this -> parent -> get_wrapped_shape_data();
	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();
	std::string primary_name = this -> primary_prop_combo_box -> currentText().toStdString();
	mass_properties -> SetInputData(shape_data[primary_name] -> get_polydata());

	mass_properties -> Update();



	std::vector<arma::vec> position_from_file, attitude_from_file;
	std::vector<double> position_time_from_file, attitude_time_from_file;


	if (shape_properties_widget -> position_from_file_button -> isChecked()){
		// This is where the text file should be parsed
		// The in the file should be
		// time (s) -- x_pos (m)-- y_pos (m) -- z_pos (m)-- x_dot (m) -- y_dot (m) -- z_dot
		ObsWindow::load_state_from_file(
			shape_properties_widget -> get_position_state_file_path(),
			position_from_file,
			position_time_from_file);
	}

	if (shape_properties_widget -> attitude_from_file_button -> isChecked()){
		// This is where the text file should be parsed
		// time (s) -- sigma_1 -- sigma_2 -- sigma_3 -- omega_x -- omega_y -- omega_z
		ObsWindow::load_state_from_file(shape_properties_widget -> get_attitude_state_file_path(),
			attitude_from_file,attitude_time_from_file);
		
	}

	for (int i  = 0; i < this -> N_images_sbox -> value(); ++i){

		double time = imaging_times[i];

		if (shape_properties_widget -> position_from_keplerian_button -> isChecked()){

			// If the considered shape is a primary, then its position state should be zero (both in position
			// and velocity)

			if (shape_properties_widget -> isPrimary()){
				positions_vec[i].push_back(arma::zeros<arma::vec>(3));
				velocities_vec[i].push_back(arma::zeros<arma::vec>(3));
			}
			else{
				arma::vec elements = shape_properties_widget -> get_orbital_elements();
				double mu = arma::datum::G * shape_properties_widget -> get_density() * mass_properties -> GetVolume();

				OC::KepState kep_state(elements,mu);
				OC::CartState cart_state = kep_state.convert_to_cart(time);
				positions_vec[i].push_back(cart_state.get_position_vector());
				velocities_vec[i].push_back(cart_state.get_velocity_vector());
			}

		} else {

			arma::vec position_state = ObsWindow::interpolate(time,position_from_file,position_time_from_file,false);

			positions_vec[i].push_back(position_state.rows(0,2));
			velocities_vec[i].push_back(position_state.rows(3,5));

		}

		if (shape_properties_widget -> attitude_from_simple_spin_button -> isChecked()){
			double period = shape_properties_widget -> get_period();
			double w = 2 * arma::datum::pi / period;
			double phi = w * time;
			phi = phi - int(phi / (2 * arma::datum::pi) ) * (2 * arma::datum::pi);

			arma::vec spin = shape_properties_widget -> get_spin();
			arma::vec mrp;

			if (phi > arma::datum::pi){
				phi = 2 * arma::datum::pi - phi;
				spin = -spin;
			}

			mrp = std::tan(phi / 4) * spin;
			arma::vec omega = w * spin;
			mrps_vec[i].push_back(mrp);
			omegas_vec[i].push_back(omega);

		} else{

			arma::vec attitude_state = ObsWindow::interpolate(time,attitude_from_file,attitude_time_from_file,true);

			mrps_vec[i].push_back(attitude_state.rows(0,2));
			omegas_vec[i].push_back(attitude_state.rows(3,5));

		}

	}

}


void ObsWindow::load_state_from_file(const std::string & filepath,
	std::vector<arma::vec> & state_vec,
	std::vector<double> & time_vec){

	std::ifstream file(filepath);
	double time, x,y,z,X,Y,Z;

	while(file >> time >> x >> y >> z >> X >> Y >> Z ){
		time_vec.push_back(time);
		arma::vec state = {x,y,z,X,Y,Z};
		state_vec.push_back(state);
	}

}

arma::vec ObsWindow::interpolate(const double & time,
	const std::vector<arma::vec> & state_from_file,
	const std::vector<double> & state_time_from_file,
	bool is_attitude_state){


	int time_index = -1;
	double t0,t1;
	arma::vec f0,f1;

	for (unsigned int i = 0; i + 1 < state_time_from_file.size(); ++i){

		if (state_time_from_file[i] < time && state_time_from_file[i + 1] >= time){
			time_index = i;
			t0 = state_time_from_file[i];
			t1 = state_time_from_file[i + 1];
			f0 = state_from_file[i];
			f1 = state_from_file[i+1];

			break;
		}

	}

	if (time_index < 0){
		throw (std::runtime_error("The prescribed time does not fall within the interpolation time table"));
	}

	// If we are dealing with an attitude state, must switch MRP if need be
	// switching is detected if two consecutive MRPs yield a relative principal rotation vector of more than 
	// 160 degrees in magnitude (a switching over the sigma^Tsigma = 1 surface brings exactly 180 deg of separation)
	if (is_attitude_state){

		double prv_relative_norm = arma::norm(RBK::dcm_to_prv(RBK::mrp_to_dcm(f0.rows(0,2))* RBK::mrp_to_dcm(-f1.rows(0,2))));
		
		if (prv_relative_norm * 180./arma::datum::pi > 90 ){
			f0.rows(0,2) = RBK::shadow_mrp(f0.rows(0,2),true);
		}


	}


	arma::vec interpolant = (time - t0)/(t1 - t0) * (f1 - f0) + f0;


	return interpolant;


}

