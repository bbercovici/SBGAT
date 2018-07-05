#include "ShapePropertiesWidget.hpp"

using namespace SBGAT_GUI;

ShapePropertiesWidget::ShapePropertiesWidget(ObsWindow * parent,bool is_primary,
	std::string title){

	QGridLayout * shape_properties_group_layout = new QGridLayout(this);

	this -> spin_raan_label = new QLabel("Spin RAAN (deg)",this);
	this -> spin_inc_label = new QLabel("Spin inclination (deg)",this);
	this -> period_label = new QLabel("Rotation period (hours)",this);


	this -> setTitle(QString::fromStdString(title));

	this -> spin_raan_sbox = new QDoubleSpinBox(this);
	this -> spin_inc_sbox = new QDoubleSpinBox(this);
	this -> period_sbox = new QDoubleSpinBox(this);

	shape_properties_group_layout -> addWidget(spin_raan_label,0,0,1,1);
	shape_properties_group_layout -> addWidget(this -> spin_raan_sbox,0,1,1,1);

	shape_properties_group_layout -> addWidget(spin_inc_label,1,0,1,1);
	shape_properties_group_layout -> addWidget(this -> spin_inc_sbox,1,1,1,1);

	shape_properties_group_layout -> addWidget(period_label,2,0,1,1);
	shape_properties_group_layout -> addWidget(this -> period_sbox,2,1,1,1);

	if (is_primary){
		this -> density_label = new QLabel("Density (kg/m^3)",this);
		this -> density_sbox = new QDoubleSpinBox(this);
		shape_properties_group_layout -> addWidget(this -> density_label,3,0,1,1);
		shape_properties_group_layout -> addWidget(this -> density_sbox,3,1,1,1);
	}

	this -> init();
	
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

	if (this -> density_sbox != nullptr){
		this -> density_sbox -> setRange(1e-10,1e10);
		this -> density_sbox -> setDecimals(6);
		this -> density_sbox -> setValue(2000);
	}
	

	this -> period_sbox -> setValue(1);


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

