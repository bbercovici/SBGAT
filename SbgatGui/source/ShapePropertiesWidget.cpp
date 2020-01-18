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



#include "ShapePropertiesWidget.hpp"
#include "SurfacePGMWindow.hpp"
#include "MassPropertiesWindow.hpp"
#include <QRadioButton>

using namespace SBGAT_GUI;

ShapePropertiesWidget::ShapePropertiesWidget(ObsWindow * parent,bool is_primary,std::string title){
	
	this -> is_primary = is_primary;
	QGridLayout * shape_properties_group_layout = new QGridLayout(this);


	this -> position_box = new QGroupBox("Position and velocity",this);
	this -> attitude_box = new QGroupBox("Attitude and angular velocity",this);


	this -> position_keplerian_widget = new QWidget(this -> position_box);
	this -> attitude_simple_spin_widget = new QWidget(this -> attitude_box);

	this -> position_from_file_widget = new QWidget(this -> position_box);
	this -> attitude_from_file_widget = new QWidget(this -> attitude_box);

	this -> load_position_state_from_file_button = new QPushButton("Select file",this -> position_from_file_widget);
	this -> load_attitude_state_from_file_button = new QPushButton("Select file",this -> attitude_from_file_widget);

	QWidget * button_box_widget_position = new QWidget(this -> position_box);
	QHBoxLayout * button_box_widget_position_layout = new QHBoxLayout(button_box_widget_position);

	QWidget * button_box_widget_attitude = new QWidget(this -> attitude_box);
	QHBoxLayout * button_box_widget_attitude_layout = new QHBoxLayout(button_box_widget_attitude);

	this -> position_from_keplerian_button = new QRadioButton(tr("Keplerian"));
	this -> position_from_file_button = new QRadioButton(tr("From file"));
	
	button_box_widget_position_layout -> addWidget(this -> position_from_keplerian_button);
	button_box_widget_position_layout -> addWidget(this -> position_from_file_button);

	this -> attitude_from_simple_spin_button = new QRadioButton(tr("Simple spin"));
	this -> attitude_from_file_button = new QRadioButton(tr("From file"));
	
	button_box_widget_attitude_layout -> addWidget(this -> attitude_from_simple_spin_button);
	button_box_widget_attitude_layout -> addWidget(this -> attitude_from_file_button);

	QGridLayout * position_box_layout = new QGridLayout(this -> position_box);
	QGridLayout * attitude_box_layout = new QGridLayout(this -> attitude_box);


	QGridLayout * position_keplerian_widget_layout = new QGridLayout(this -> position_keplerian_widget);
	QGridLayout * attitude_simple_spin_widget_layout = new QGridLayout(this -> attitude_simple_spin_widget);

	QGridLayout * position_from_file_widget_layout = new QGridLayout(this -> position_from_file_widget);
	QGridLayout * attitude_from_file_widget_layout = new QGridLayout(this -> attitude_from_file_widget);


	this -> load_position_state_from_file_label = new QLabel("(No file selected)",this -> position_from_file_widget);
	this -> load_attitude_state_from_file_label = new QLabel("(No file selected)",this -> attitude_from_file_widget);

	position_from_file_widget_layout -> addWidget(this -> load_position_state_from_file_button,0,0,1,1);
	position_from_file_widget_layout -> addWidget(this -> load_position_state_from_file_label,0,1,1,1);


	attitude_from_file_widget_layout -> addWidget(this -> load_attitude_state_from_file_button,0,0,1,1);
	attitude_from_file_widget_layout -> addWidget(this -> load_attitude_state_from_file_label,0,1,1,1);


	position_box_layout -> addWidget(button_box_widget_position);
	position_box_layout -> addWidget(this -> position_keplerian_widget);
	position_box_layout -> addWidget(this -> position_from_file_widget);


	attitude_box_layout -> addWidget(button_box_widget_attitude);
	attitude_box_layout -> addWidget(this -> attitude_simple_spin_widget);
	attitude_box_layout -> addWidget(this -> attitude_from_file_widget);

	shape_properties_group_layout -> addWidget(this -> attitude_box,1,0,1,2);


	this -> spin_raan_label = new QLabel("Spin RAAN (deg)",this);
	this -> spin_inc_label = new QLabel("Spin inclination (deg)",this);
	this -> period_label = new QLabel("Rotation period (hours)",this);



	this -> setTitle(QString::fromStdString(title));

	this -> spin_raan_sbox = new QDoubleSpinBox(this);
	this -> spin_inc_sbox = new QDoubleSpinBox(this);
	this -> period_sbox = new QDoubleSpinBox(this);

	attitude_simple_spin_widget_layout -> addWidget(spin_raan_label,0,0,1,1);
	attitude_simple_spin_widget_layout -> addWidget(this -> spin_raan_sbox,0,1,1,1);

	attitude_simple_spin_widget_layout -> addWidget(spin_inc_label,1,0,1,1);
	attitude_simple_spin_widget_layout -> addWidget(this -> spin_inc_sbox,1,1,1,1);

	attitude_simple_spin_widget_layout -> addWidget(period_label,2,0,1,1);
	attitude_simple_spin_widget_layout -> addWidget(this -> period_sbox,2,1,1,1);


	

	shape_properties_group_layout -> addWidget(this -> position_box,0,0,1,2);

	if (is_primary){
		this -> position_from_keplerian_button -> setText("Fixed");
	}
	else{
		QLabel * sma_label = new QLabel("Semi-major axis (m)",this);
		QLabel * ecc_label = new QLabel("Eccentricity",this);
		QLabel * inc_label = new QLabel("Inclination (deg)",this);
		QLabel * Omega_label = new QLabel("RAAN (deg)",this);
		QLabel * omega_label = new QLabel("Longitude of perigee (deg)",this);
		QLabel * M0_label = new QLabel("Mean anomaly at epoch (deg)",this);
		QLabel * density_label = new QLabel("Density of central body (kg/m^3)",this);

		this -> sma_sbox  = new QDoubleSpinBox(this);
		this -> ecc_sbox  = new QDoubleSpinBox(this);
		this -> inc_sbox  = new QDoubleSpinBox(this);
		this -> Omega_sbox  = new QDoubleSpinBox(this);
		this -> omega_sbox  = new QDoubleSpinBox(this);
		this -> M0_sbox  = new QDoubleSpinBox(this);
		this -> density_sbox = new QDoubleSpinBox(this);

		position_keplerian_widget_layout -> addWidget(sma_label,0,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> sma_sbox,0,1,1,1);

		position_keplerian_widget_layout -> addWidget(ecc_label,1,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> ecc_sbox,1,1,1,1);

		position_keplerian_widget_layout -> addWidget(inc_label,2,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> inc_sbox,2,1,1,1);

		position_keplerian_widget_layout -> addWidget(Omega_label,3,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> Omega_sbox,3,1,1,1);

		position_keplerian_widget_layout -> addWidget(omega_label,4,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> omega_sbox,4,1,1,1);

		position_keplerian_widget_layout -> addWidget(M0_label,5,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> M0_sbox,5,1,1,1);

		position_keplerian_widget_layout -> addWidget(density_label,6,0,1,1);
		position_keplerian_widget_layout -> addWidget(this -> density_sbox,6,1,1,1);

	}

	this -> init();

	connect(this -> position_from_keplerian_button, SIGNAL(clicked()),this,SLOT(toggle_position_visibility()));
	connect(this -> position_from_file_button, SIGNAL(clicked()),this,SLOT(toggle_position_visibility()));

	connect(this -> attitude_from_simple_spin_button, SIGNAL(clicked()),this,SLOT(toggle_attitude_visibility()));
	connect(this -> attitude_from_file_button, SIGNAL(clicked()),this,SLOT(toggle_attitude_visibility()));

	connect(this -> load_position_state_from_file_button,SIGNAL(clicked()),this,SLOT(load_position_state()));
	connect(this -> load_attitude_state_from_file_button,SIGNAL(clicked()),this,SLOT(load_attitude_state()));

	
}



ShapePropertiesWidget::ShapePropertiesWidget(SurfacePGMWindow * parent,std::string title){
	
	QGridLayout * attitude_box_layout = new QGridLayout(this);

	this -> spin_raan_label = new QLabel("Spin RAAN (deg)",this);
	this -> spin_inc_label = new QLabel("Spin inclination (deg)",this);
	this -> period_label = new QLabel("Rotation period (hours)",this);
	QLabel * density_label = new QLabel("Density of central body (kg/m^3)",this);

	this -> setTitle(QString::fromStdString(title));

	this -> spin_raan_sbox = new QDoubleSpinBox(this);
	this -> spin_inc_sbox = new QDoubleSpinBox(this);
	this -> period_sbox = new QDoubleSpinBox(this);
	this -> density_sbox = new QDoubleSpinBox(this);

	attitude_box_layout -> addWidget(spin_raan_label,0,0,1,1);
	attitude_box_layout -> addWidget(this -> spin_raan_sbox,0,1,1,1);

	attitude_box_layout -> addWidget(spin_inc_label,1,0,1,1);
	attitude_box_layout -> addWidget(this -> spin_inc_sbox,1,1,1,1);

	attitude_box_layout -> addWidget(period_label,2,0,1,1);
	attitude_box_layout -> addWidget(this -> period_sbox,2,1,1,1);
	
	attitude_box_layout -> addWidget(density_label,3,0,1,1);
	attitude_box_layout -> addWidget(this -> density_sbox,3,1,1,1);

	this -> spin_raan_sbox -> setDecimals(6);
	this -> spin_inc_sbox -> setDecimals(6);

	this -> period_sbox -> setDecimals(6);

	this -> spin_raan_sbox -> setRange(-180,180);
	this -> spin_inc_sbox -> setRange(-180,180);
	this -> period_sbox -> setRange(1e-10,1e10);
	this -> period_sbox -> setValue(12);


	this -> spin_raan_sbox -> setValue(0);
	this -> spin_inc_sbox -> setValue(0);

	this -> density_sbox -> setRange(1e-10,1e10);
	this -> density_sbox -> setDecimals(6);
	this -> density_sbox -> setValue(2000);
}


ShapePropertiesWidget::ShapePropertiesWidget(MassPropertiesWindow * parent,std::string title){
	
	QGridLayout * attitude_box_layout = new QGridLayout(this);

	QLabel * density_label = new QLabel("Density of central body (kg/m^3)",this);

	this -> setTitle(QString::fromStdString(title));

	this -> density_sbox = new QDoubleSpinBox(this);
	
	attitude_box_layout -> addWidget(density_label,0,0,1,1);
	attitude_box_layout -> addWidget(this -> density_sbox,0,1,1,1);


	this -> density_sbox -> setRange(1e-10,1e10);
	this -> density_sbox -> setDecimals(6);
	this -> density_sbox -> setValue(2000);
}





void ShapePropertiesWidget::init(){

	this -> spin_raan_sbox -> setDecimals(6);
	this -> spin_inc_sbox -> setDecimals(6);

	this -> period_sbox -> setDecimals(6);

	this -> spin_raan_sbox -> setRange(-180,180);
	this -> spin_inc_sbox -> setRange(-180,180);
	this -> period_sbox -> setRange(1e-10,1e10);

	this -> spin_raan_sbox -> setValue(0);
	this -> spin_inc_sbox -> setValue(0);

	if(!this -> is_primary){
		this -> sma_sbox -> setRange(1e-10,1e10);
		this -> ecc_sbox -> setRange(0,1);
		this -> inc_sbox -> setRange(0,180);
		this -> Omega_sbox -> setRange(-180,180);
		this -> omega_sbox -> setRange(-180,180);
		this -> M0_sbox -> setRange(-180,180);


		this -> sma_sbox -> setValue(1000);
		this -> ecc_sbox -> setValue(0);
		this -> inc_sbox -> setValue(0);
		this -> Omega_sbox -> setValue(0);
		this -> omega_sbox -> setValue(0);
		this -> M0_sbox -> setValue(0);

		this -> sma_sbox -> setDecimals(6);
		this -> ecc_sbox -> setDecimals(6);
		this -> inc_sbox -> setDecimals(6);
		this -> Omega_sbox -> setDecimals(6);
		this -> omega_sbox -> setDecimals(6);
		this -> M0_sbox -> setDecimals(6);

		this -> density_sbox -> setRange(1e-10,1e10);
		this -> density_sbox -> setDecimals(6);
		this -> density_sbox -> setValue(2000);

	}

	this -> position_from_keplerian_button -> setChecked(1);
	this -> attitude_from_simple_spin_button -> setChecked(1);

	this -> period_sbox -> setValue(1);
	this -> position_from_file_widget -> setVisible(0);
	this -> attitude_from_file_widget -> setVisible(0);
	
}


arma::vec ShapePropertiesWidget::get_spin() const{

	double d2r = arma::datum::pi /180;
	arma::vec spin = {0,0,1};
	spin = (RBK::M2(this -> spin_inc_sbox -> value() * d2r) 
		* RBK::M3(this -> spin_raan_sbox -> value() * d2r)).t() * spin;
	return spin;

}

double ShapePropertiesWidget::get_period() const{
	return this -> period_sbox -> value() * 3600; 
}

double ShapePropertiesWidget::get_density() const{

	if(this -> density_sbox){
		return this -> density_sbox -> value();
	}
	else{
		return 0;
	}

}


void ShapePropertiesWidget::toggle_position_visibility(){

	if (this -> position_from_keplerian_button -> isChecked()){
		
		this -> position_keplerian_widget -> setVisible(1);
		this -> position_from_file_widget -> setVisible(0);

		this -> repaint();
		

	}
	else{
		this -> position_keplerian_widget -> setVisible(0);
		this -> position_from_file_widget -> setVisible(1);

		this -> repaint();
	}

}

void ShapePropertiesWidget::toggle_attitude_visibility(){

	if (this -> attitude_from_simple_spin_button -> isChecked()){
		this -> attitude_simple_spin_widget -> setVisible(1);
		this -> attitude_from_file_widget -> setVisible(0);

		this -> repaint();

	}
	else{
		this -> attitude_simple_spin_widget -> setVisible(0);
		this -> attitude_from_file_widget -> setVisible(1);

		this -> repaint();
	}

}

arma::vec ShapePropertiesWidget::get_orbital_elements() const{
	double d2r = arma::datum::pi /180;

	arma::vec orbital_elements(6);

	orbital_elements(0) = this -> sma_sbox -> value();
	orbital_elements(1) = this -> ecc_sbox -> value();
	orbital_elements(2) = this -> inc_sbox -> value() * d2r;
	orbital_elements(3) = this -> Omega_sbox -> value() * d2r;
	orbital_elements(4) = this -> omega_sbox -> value() * d2r;
	orbital_elements(5) = this -> M0_sbox-> value() * d2r;

	return orbital_elements;

}



void ShapePropertiesWidget::load_position_state(){

	QString fileName = QFileDialog::getOpenFileName(this, tr("Load trajectory file"),
		"/home",
		tr(" (*.txt)"));
	
	if (!fileName.isNull()){
		this -> load_position_state_from_file_label -> setText(fileName);
		this -> position_state_file_path = fileName.toStdString();
		
	}

}
void ShapePropertiesWidget::load_attitude_state(){

	QString fileName = QFileDialog::getOpenFileName(this, tr("Load attitude history file"),
		"/home",
		tr(" (*.txt)"));

	if (!fileName.isNull()){
		this -> load_attitude_state_from_file_label -> setText(fileName);
		this -> attitude_state_file_path = fileName.toStdString();
	}




}


std::string ShapePropertiesWidget::get_position_state_file_path() const{
	return this -> position_state_file_path;
}

std::string ShapePropertiesWidget::get_attitude_state_file_path() const{
	return this -> attitude_state_file_path;
}







