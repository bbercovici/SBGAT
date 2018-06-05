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

#include "YORPWindow.hpp"
#include <SBGATSrpYorp.hpp>

using namespace SBGAT_GUI;

YORPWindow::YORPWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("YORP Coefficients Computation");


	QVBoxLayout * yorp_layout = new QVBoxLayout(this);

	QGroupBox * shape_model_group = new QGroupBox(tr("Ray-traced shape"));
	QGroupBox * surface_group = new QGroupBox(tr("Surface properties"));
	QGroupBox * fidelity_group = new QGroupBox(tr("Fidelity"));
	QGroupBox * fourier_group = new QGroupBox(tr("Fourier expansion"));
	QGroupBox * output_group = new QGroupBox(tr("Output directory"));


	
	QGridLayout * shape_model_group_layout = new QGridLayout(shape_model_group);
	QGridLayout * surface_group_layout = new QGridLayout(surface_group);
	QGridLayout * fidelity_group_layout = new QGridLayout(fidelity_group);
	QGridLayout * fourier_group_layout = new QGridLayout(fourier_group);
	QGridLayout * output_group_layout = new QGridLayout(output_group);



	QLabel * shape_label = new QLabel("Shape model",this);
	QLabel * rho_label = new QLabel("Rho",this);
	QLabel * spec_label = new QLabel("Specular",this);
	QLabel * bounces_label = new QLabel("Reflexions",this);
	QLabel * refine_label = new QLabel("Refinement",this);
	QLabel * order_label = new QLabel("Maximum order",this);
	QLabel * voxel_label = new QLabel("Voxels per dimension",this);


	QLabel * lambda_del_label = new QLabel("Azimuth step (degree)",this);
	QLabel * delta_del_label = new QLabel("Declination step (degree)",this);


	this -> prop_combo_box = new QComboBox (this);

	this -> fourier_combo_box = new QComboBox (this);
	this -> bounces_combo_box = new QComboBox (this);
	this -> refine_combo_box = new QComboBox (this);
	this -> voxel_combo_box = new QComboBox (this);

	this -> rho_sbox  = new QDoubleSpinBox(this);
	this -> spec_sbox  = new QDoubleSpinBox(this);

	this -> lambda_del_sbox = new QDoubleSpinBox(this);
	this -> delta_del_sbox = new QDoubleSpinBox(this);

	this -> rho_sbox -> setRange(0,1);
	this -> spec_sbox -> setRange(0,1);
	
	this -> lambda_del_sbox -> setRange(1e-3,10);
	this -> delta_del_sbox -> setRange(1e-3,10);

	this ->  open_output_file_dialog_button = new QPushButton("Select output directory",this);

	shape_model_group_layout -> addWidget(shape_label,0,0,1,1);
	shape_model_group_layout -> addWidget(this -> prop_combo_box,0,1,1,1);

	surface_group_layout -> addWidget(rho_label,0,0,1,1);
	surface_group_layout -> addWidget(this -> rho_sbox,0,1,1,1);

	surface_group_layout -> addWidget(spec_label,1,0,1,1);
	surface_group_layout -> addWidget(this -> spec_sbox,1,1,1,1);

	fidelity_group_layout -> addWidget(bounces_label,0,0,1,1);
	fidelity_group_layout -> addWidget(this -> bounces_combo_box,0,1,1,1);
	fidelity_group_layout -> addWidget(refine_label,1,0,1,1);
	fidelity_group_layout -> addWidget(this -> refine_combo_box,1,1,1,1);
	fidelity_group_layout -> addWidget(voxel_label,2,0,1,1);
	fidelity_group_layout -> addWidget(this -> voxel_combo_box,2,1,1,1);

	fidelity_group_layout -> addWidget(lambda_del_label,3,0,1,1);
	fidelity_group_layout -> addWidget(this -> lambda_del_sbox,3,1,1,1);


	fidelity_group_layout -> addWidget(delta_del_label,4,0,1,1);
	fidelity_group_layout -> addWidget(this -> delta_del_sbox,4,1,1,1);

	fourier_group_layout -> addWidget(order_label,0,0,1,1);
	fourier_group_layout -> addWidget(this -> fourier_combo_box,0,1,1,1);

	output_group_layout -> addWidget(this -> open_output_file_dialog_button, 0, 0, 1, 2);


	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);

	yorp_layout -> addWidget(shape_model_group);
	yorp_layout -> addWidget(surface_group);
	yorp_layout -> addWidget(fidelity_group);
	yorp_layout -> addWidget(fourier_group);
	yorp_layout -> addWidget(output_group);
	yorp_layout -> addWidget(button_box);

	this -> init();


	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> open_output_file_dialog_button,SIGNAL(clicked()),this,
		SLOT(open_output_file_dialog()));

}

void YORPWindow::init(){

	this -> rho_sbox -> setValue(1);
	this -> lambda_del_sbox -> setValue(1.);
	this -> delta_del_sbox -> setValue(1.);

	for (unsigned int i = 0; i < 10 ; ++i) {
		this -> fourier_combo_box -> insertItem(i, QString::number(i));
	}

	this -> fourier_combo_box -> setCurrentIndex(2);

	for (unsigned int i = 0; i < 10 ; ++i) {
		this -> bounces_combo_box -> insertItem(i, QString::number(i));
	}

	this -> bounces_combo_box -> setCurrentIndex(3);

	for (unsigned int i = 0; i < 10 ; ++i) {
		this -> refine_combo_box -> insertItem(i, QString::number(i));
	}
	
	this -> refine_combo_box -> setCurrentIndex(5);


	for (unsigned int i = 0; i < 15 ; ++i) {
		this -> voxel_combo_box -> insertItem(i, QString::number(5 * (i + 1) ));
	}

	this -> voxel_combo_box -> setCurrentIndex(7);



	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();
	
	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
	}
	
	this -> button_box -> button(QDialogButtonBox::Ok) -> setEnabled(false);





}

void YORPWindow::open_output_file_dialog(){

	QString path = QFileDialog::getExistingDirectory(this, tr("Save directory"), "~");
	this -> output_path = path.toStdString();
	if (this -> output_path.size() > 0){
		this -> button_box -> button(QDialogButtonBox::Ok) -> setEnabled(true);

	}
}


void YORPWindow::accept(){

	vtkSmartPointer<SBGATSrpYorp> yorp = vtkSmartPointer<SBGATSrpYorp>::New();

	std::string name = this -> prop_combo_box -> currentText().toStdString();

	auto shape_data = this -> parent -> get_wrapped_shape_data();

	yorp -> SetInputData(shape_data[name] -> get_polydata());

	

	yorp -> set_rho(this -> rho_sbox -> value());
	yorp -> set_spec(this -> spec_sbox -> value());
	yorp -> set_lambdaDel(this -> lambda_del_sbox -> value());
	yorp -> set_deltaDel(this -> delta_del_sbox -> value());
	yorp -> set_howManyBounces(this -> bounces_combo_box -> currentText().toInt());
	yorp -> set_maxFourier(this -> fourier_combo_box -> currentText().toInt());
	yorp -> set_numrefine(this -> refine_combo_box -> currentText().toInt());
	yorp -> set_numVox(this -> voxel_combo_box -> currentText().toInt());


	if (this -> output_path.size() == 0){
		yorp -> set_outputFileBaseName("./output.txt");
	}
	else{
		yorp -> set_outputFileBaseName(this -> output_path);
	}

	yorp -> Update();




	QDialog::accept();
}


