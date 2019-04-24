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
#include <QHeaderView>

using namespace SBGAT_GUI;

ShapeUncertaintyWidget::ShapeUncertaintyWidget(SurfacePGMWindow * parent,std::string title){
	
	this -> parent = parent;
	QGridLayout * shape_uncertainty_group_layout = new QGridLayout(this);
	this -> setTitle(QString::fromStdString(title));

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
	QLabel * global_covariance_regularization_label = new QLabel("Number of covariance regularizations");
	QLabel * local_covariance_regularization_label = new QLabel("Number of covariance regularizations");


	this -> add_region_button = new QPushButton("Add Region",this);
	this -> remove_region_button = new QPushButton("Remove Region",this);

	this -> global_sigma_spinbox = new QDoubleSpinBox(this);
	this -> global_correlation_distance_spinbox = new QDoubleSpinBox(this);

	this -> global_sigma_spinbox -> setMinimum(0);
	this -> global_correlation_distance_spinbox -> setMinimum(1e-2);

	this -> global_sigma_spinbox -> setValue(0);
	this -> global_correlation_distance_spinbox -> setValue(0);


	this -> global_covariance_regularization_spin_box = new QSpinBox(uncertainty_global_widget);
	this -> global_covariance_regularization_spin_box -> setMinimum(0);
	this -> global_covariance_regularization_spin_box -> setValue(2);

	this -> local_covariance_regularization_spin_box = new QSpinBox(uncertainty_local_widget);
	this -> local_covariance_regularization_spin_box -> setMinimum(0);
	this -> local_covariance_regularization_spin_box -> setValue(2);
	
	this -> local_uncertainty_table_widget = new QTableWidget(this);

	this -> local_uncertainty_table_widget -> setColumnCount(3);

	QTableWidgetItem * center_facet_index_header = new QTableWidgetItem("Center Facet Index");
	QTableWidgetItem * sigma_header = new QTableWidgetItem("Noise standard-deviation (m)");
	QTableWidgetItem * correlation_length_header = new QTableWidgetItem("Correlation length(m)");

	this -> local_uncertainty_table_widget -> setHorizontalHeaderItem(0, center_facet_index_header);
	this -> local_uncertainty_table_widget -> setHorizontalHeaderItem(1, sigma_header);
	this -> local_uncertainty_table_widget -> setHorizontalHeaderItem(2, correlation_length_header);
	this -> local_uncertainty_table_widget -> horizontalHeader()->setSectionResizeMode(QHeaderView::Stretch);


	uncertainty_global_layout -> addWidget(global_sigma_label,0,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_sigma_spinbox,0,1,1,1);

	uncertainty_global_layout -> addWidget(global_correlation_distance_label,1,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_correlation_distance_spinbox,1,1,1,1);

	uncertainty_global_layout -> addWidget(global_covariance_regularization_label,2,0,1,1);
	uncertainty_global_layout -> addWidget(this -> global_covariance_regularization_spin_box,2,1,1,1);

	uncertainty_local_layout -> addWidget(this -> local_uncertainty_table_widget,0,0,3,6);
	uncertainty_local_layout -> addWidget(this -> add_region_button,3,0,1,3);
	uncertainty_local_layout -> addWidget(this -> remove_region_button,3,3,1,3);
	uncertainty_local_layout -> addWidget(local_covariance_regularization_label,4,0,1,3);
	uncertainty_local_layout -> addWidget(this -> local_covariance_regularization_spin_box,4,3,1,3);

	connect(this -> covariance_input_file_button,SIGNAL(clicked()),this,SLOT(select_covariance_input_file()));
	connect(this -> add_region_button,SIGNAL(clicked()),this,SLOT(add_shape_uncertainty_region()));
	connect(this -> remove_region_button,SIGNAL(clicked()),this,SLOT(remove_shape_uncertainty_region()));

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

void ShapeUncertaintyWidget::add_shape_uncertainty_region(){
	QSpinBox * facet_center_index_spinbox = new QSpinBox(this -> local_uncertainty_table_widget);
	facet_center_index_spinbox -> setMinimum(0);
	facet_center_index_spinbox -> setMaximum(static_cast<int>(this -> parent -> parent -> get_wrapped_shape_data()[this -> parent -> primary_prop_combo_box -> currentText().toStdString()] -> get_polydata() -> GetNumberOfCells()) - 1);

	QDoubleSpinBox * standard_deviation_spin_box = new QDoubleSpinBox( this -> local_uncertainty_table_widget );
	QDoubleSpinBox * correlation_distance_spin_box = new QDoubleSpinBox( this -> local_uncertainty_table_widget );

	standard_deviation_spin_box -> setMinimum(0);
	correlation_distance_spin_box -> setMinimum(0);


	this -> local_uncertainty_table_widget -> setRowCount( this -> local_uncertainty_table_widget -> rowCount() + 1 );
	this -> local_uncertainty_table_widget -> setCellWidget ( this -> local_uncertainty_table_widget -> rowCount() - 1, 0, facet_center_index_spinbox );
	this -> local_uncertainty_table_widget -> setCellWidget ( this -> local_uncertainty_table_widget -> rowCount() - 1, 1, standard_deviation_spin_box );
	this -> local_uncertainty_table_widget -> setCellWidget ( this -> local_uncertainty_table_widget -> rowCount() - 1, 2, correlation_distance_spin_box );



}

void ShapeUncertaintyWidget::remove_shape_uncertainty_region(){

	if (this -> local_uncertainty_table_widget -> rowCount() > 0){
		this -> local_uncertainty_table_widget -> setRowCount( this -> local_uncertainty_table_widget -> rowCount() - 1 );
	}
}

int ShapeUncertaintyWidget::get_global_covariance_regularization_number() const{

	return this -> global_covariance_regularization_spin_box -> value();
}

int ShapeUncertaintyWidget::get_local_covariance_regularization_number() const{

	return this -> local_covariance_regularization_spin_box -> value();
}	

void ShapeUncertaintyWidget::clear(){
	this -> local_uncertainty_table_widget -> setRowCount(0);
}




