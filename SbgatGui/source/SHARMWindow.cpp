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

#include "SHARMWindow.hpp"
#include <SBGATSphericalHarmo.hpp>

using namespace SBGAT_GUI;

SHARMWindow::SHARMWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Gravity Spherical Harmonics Computation");


	QVBoxLayout * sharm_layout = new QVBoxLayout(this);

	QGroupBox * shape_model_group = new QGroupBox(tr("Shape"));
	QGroupBox * properties_group = new QGroupBox(tr("Shape properties"));
	QGroupBox * settings_group = new QGroupBox(tr("Expansion settings"));
	QGroupBox * output_group = new QGroupBox(tr("Output directory"));

	
	QGridLayout * shape_model_group_layout = new QGridLayout(shape_model_group);
	QGridLayout * settings_group_layout = new QGridLayout(settings_group);
	QGridLayout * properties_group_layout = new QGridLayout(properties_group);
	QGridLayout * output_group_layout = new QGridLayout(output_group);

	QWidget * choose_normalization_widget = new QWidget(this);


	this -> normalized_button = new QRadioButton("Normalized", choose_normalization_widget);
	this -> non_normalized_button = new QRadioButton("Non-Normalized", choose_normalization_widget);

	QLabel * shape_label = new QLabel("Shape model",this);
	QLabel * density_label = new QLabel("Density (kg/m^3)",this);
	QLabel * ref_radius_label = new QLabel("Reference radius (m)",this);
	QLabel * degree_label = new QLabel("Expansion degree",this);

	this -> prop_combo_box = new QComboBox (this);

	this -> density_sbox = new QDoubleSpinBox(this);
	this -> ref_radius_sbox = new QDoubleSpinBox(this);
	this -> degree_combo_box = new QComboBox (this);

	this ->  open_output_file_dialog_button = new QPushButton("Select output directory",this);

	shape_model_group_layout -> addWidget(shape_label,0,0,1,1);
	shape_model_group_layout -> addWidget(this -> prop_combo_box,0,1,1,1);

	properties_group_layout -> addWidget(density_label,0,0,1,1);
	properties_group_layout -> addWidget(this -> density_sbox,0,1,1,1);

	properties_group_layout -> addWidget(ref_radius_label,1,0,1,1);
	properties_group_layout -> addWidget(this -> ref_radius_sbox,1,1,1,1);

	settings_group_layout -> addWidget(degree_label,0,0,1,1);
	settings_group_layout -> addWidget(this -> degree_combo_box,0,1,1,1);
	settings_group_layout -> addWidget(this -> normalized_button,1,0,1,1);
	settings_group_layout -> addWidget(this -> non_normalized_button,1,1,1,1);

	output_group_layout -> addWidget(this -> open_output_file_dialog_button, 0, 0, 1, 2);


	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);

	sharm_layout -> addWidget(shape_model_group);
	sharm_layout -> addWidget(properties_group);
	sharm_layout -> addWidget(settings_group);
	sharm_layout -> addWidget(output_group);
	sharm_layout -> addWidget(button_box);


	connect(this -> prop_combo_box, SIGNAL(currentIndexChanged(int)),this,SLOT(changed_selected_shape()));
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> open_output_file_dialog_button,SIGNAL(clicked()),this,
		SLOT(open_output_file_dialog()));

	this -> init();


}

void SHARMWindow::init(){

	this -> density_sbox -> setRange(0,1e10);
	this -> ref_radius_sbox -> setRange(0,1e10);


	this -> density_sbox -> setValue(2000.);
	this -> ref_radius_sbox -> setValue(1.);

	for (unsigned int i = 0; i < 40 ; ++i) {
		this -> degree_combo_box -> insertItem(i, QString::number(i));
	}

	this -> degree_combo_box -> setCurrentIndex(0);


	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();
	

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
	}
	
	this -> button_box -> button(QDialogButtonBox::Ok) -> setEnabled(false);

	if (wrapped_shape_data.size() == 0){
		this -> open_output_file_dialog_button -> setEnabled(false);
	}

	this -> normalized_button -> setChecked(true);
	this -> non_normalized_button -> setChecked(false);


}

void SHARMWindow::open_output_file_dialog(){

	QString path = QFileDialog::getExistingDirectory(this, tr("Save directory"), "~");
	this -> output_path = path.toStdString();
	if (this -> output_path.size() > 0){
		this -> button_box -> button(QDialogButtonBox::Ok) -> setEnabled(true);

	}
}

void SHARMWindow::accept(){

	vtkSmartPointer<SBGATSphericalHarmo> spherical_harmonics = vtkSmartPointer<SBGATSphericalHarmo>::New();
	std::string name = this -> prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();

	if ( shape_data.find(name)!= shape_data.end()){
		spherical_harmonics -> SetInputData(shape_data[name] -> get_polydata());
	}


	// The shape should be barycentered/principal axis aligned somewhere here


	spherical_harmonics -> SetDensity(this -> density_sbox -> value());
	spherical_harmonics -> SetScaleMeters();
	spherical_harmonics -> SetReferenceRadius(this -> ref_radius_sbox -> value());

	if (this -> normalized_button -> isChecked()){
		spherical_harmonics -> IsNormalized(); 
	}
	else{
		spherical_harmonics -> IsNonNormalized(); 
	}

	spherical_harmonics -> SetDegree(this -> degree_combo_box -> currentText().toInt());
	spherical_harmonics -> Update();

	auto Cnm = spherical_harmonics -> GetCnm();
	auto Snm = spherical_harmonics -> GetSnm();

	if (this -> output_path.size() != 0){
		Cnm.save(this -> output_path + "/Cnm.txt",arma::raw_ascii);
		Snm.save(this -> output_path + "/Snm.txt",arma::raw_ascii);
	}
	else{
		Cnm.save("Cnm.txt",arma::raw_ascii);
		Snm.save("Snm.txt",arma::raw_ascii);
	}

	QDialog::accept();
}

void SHARMWindow::changed_selected_shape(){

	std::string name = this -> prop_combo_box -> currentText().toStdString();

	auto shape_data = this -> parent -> get_wrapped_shape_data();

	vtkPolyData * polydata = shape_data[name] -> get_polydata();
	polydata -> ComputeBounds();
	double * bounds = polydata -> GetBounds();
	
	double x_max = std::max(std::abs(bounds[0]),std::abs(bounds[1]));
	double y_max = std::max(std::abs(bounds[2]),std::abs(bounds[3]));
	double z_max = std::max(std::abs(bounds[4]),std::abs(bounds[5]));

	double r_max[3] = {x_max,y_max,z_max};
	double radius_max = vtkMath::Norm(r_max);
	this -> ref_radius_sbox -> setValue(radius_max);


}

