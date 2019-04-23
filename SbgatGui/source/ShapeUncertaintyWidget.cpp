/** MIT License

Copyright (c) 2019 Benjamin Bercovici and Jay McMahon

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


#include "ShapeUncertaintyWidget.hpp"
#include "SurfacePGMWindow.hpp"
#include <QRadioButton>
#include <QButtonGroup>
#include <QDoubleSpinBox>
#include <QSpinBox>

using namespace SBGAT_GUI;

ShapeUncertaintyWidget::ShapeUncertaintyWidget(QDialog * parent,std::string title){
	
	QGridLayout * shape_uncertainty_group_layout = new QGridLayout(this);
	this -> setTitle(QString::fromStdString(title));
	QWidget * uncertainty_type_button_group_widget = new QWidget(this);

	this -> tab_widget = new QTabWidget(this);
	QWidget * no_uncertainty_widget = new QWidget(this);
	QWidget * uncertainty_from_file_widget = new QWidget(this);
	QWidget * uncertainty_global_widget = new QWidget(this);
	QWidget * uncertainty_local_widget = new QWidget(this);

	QGridLayout * uncertainty_from_file_layout = new QGridLayout(uncertainty_from_file_widget);
	QGridLayout * uncertainty_global_layout = new QGridLayout(uncertainty_global_widget);
	QGridLayout * uncertainty_local_layout = new QGridLayout(uncertainty_local_widget);

	this -> tab_widget -> addTab(no_uncertainty_widget,"No Uncertainty");
	this -> tab_widget -> addTab(uncertainty_from_file_widget,"From File");
	this -> tab_widget -> addTab(uncertainty_global_widget,"From Global Distribution");
	this -> tab_widget -> addTab(uncertainty_local_widget,"From Uncertainty Regions");

	this -> covariance_input_file_label = new QLabel("No File Selected",this);
	this -> covariance_input_file_button = new QPushButton("Select Input File",this);

	uncertainty_from_file_layout -> addWidget(this -> covariance_input_file_button,0,0,1,1);
	uncertainty_from_file_layout -> addWidget(this -> covariance_input_file_label,0,1,1,1);

	shape_uncertainty_group_layout -> addWidget(this -> tab_widget,0,0,3,6);
	QHBoxLayout * no_uncertainty_layout = new QHBoxLayout(no_uncertainty_widget);
	QLabel * no_uncertainty_label = new QLabel("The shape is treated as deterministic");
	no_uncertainty_layout -> addWidget(no_uncertainty_label);

	QLabel * global_sigma_label = new QLabel("Noise standard deviation (m)");
	QLabel * global_correlation_distance_label = new QLabel("Correlation distance (m)");
	QLabel * covariance_regularization_label = new QLabel("Number of covariance regularizations");


	this -> global_sigma_spinbox = new QDoubleSpinBox(this);
	this -> global_correlation_distance_spinbox = new QDoubleSpinBox(this);

	this -> global_sigma_spinbox -> setMinimum(0);
	this -> global_correlation_distance_spinbox -> setMinimum(1e-2);

	this -> global_sigma_spinbox -> setValue(0);
	this -> global_correlation_distance_spinbox -> setValue(0);


	this -> global_covariance_regularization_spin_box = new QSpinBox(uncertainty_global_widget);
	this -> global_covariance_regularization_spin_box -> setMinimum(0);
	this -> global_covariance_regularization_spin_box -> setValue(2);
	


	uncertainty_global_layout -> addWidget(global_sigma_label,0,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_sigma_spinbox,0,1,1,1);

	uncertainty_global_layout -> addWidget(global_correlation_distance_label,1,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_correlation_distance_spinbox,1,1,1,1);

	uncertainty_global_layout -> addWidget(covariance_regularization_label,2,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_covariance_regularization_spin_box,2,1,1,1);

	// uncertainty_local_layout -> addWidget(covariance_regularization_label,2,0,1,1);
	// uncertainty_local_layout -> addWidget(this -> global_covariance_regularization_spin_box,2,1,1,1);

	connect(this -> covariance_input_file_button,SIGNAL(clicked()),this,SLOT(select_covariance_input_file()));

}


void ShapeUncertaintyWidget::select_covariance_input_file(){

	QString fileName = QFileDialog::getOpenFileName(this,tr("Load shape covariance"), "~/", tr("Text file (*.txt)"));

	if (fileName.isEmpty() == false) {

		this -> covariance_input_file = fileName.toStdString();
		this -> covariance_input_file_label -> setText(fileName);
	}

	else{
		this -> covariance_input_file = "";
		this -> covariance_input_file_label -> setText("No File Selected");

	} 


}

int ShapeUncertaintyWidget::get_global_covariance_regularization_number() const{

	return global_covariance_regularization_spin_box -> value();
}

int ShapeUncertaintyWidget::get_local_covariance_regularization_number() const{

	throw(std::runtime_error("not implemented yet"));
}	




