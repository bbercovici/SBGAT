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

#include "SelectMapperWindow.hpp"

#include <QVBoxLayout>
#include <QHBoxLayout>

#include <QGroupBox>
#include <QPushButton>
#include <QCheckBox>

#include <QDialogButtonBox>
#include <QTableWidget>
#include <QLabel>
#include <vtkLookupTable.h>
#include <vtkFloatArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>




using namespace SBGAT_GUI;


SelectMapperWindow::SelectMapperWindow(Mainwindow * parent) : QDialog(parent) {

	this -> parent = parent;
	this -> setWindowTitle("Select Shape Mapper");

	QVBoxLayout * main_layout = new QVBoxLayout(this);

	this -> prop_combo_box = new QComboBox (this);
	this -> mapper_combo_box = new QComboBox (this);


	QWidget * shape_widget = new QWidget(this);
	QHBoxLayout * shape_widget_layout = new QHBoxLayout(shape_widget);
	shape_widget_layout -> addWidget(new QLabel("Shape",this));
	shape_widget_layout -> addWidget(this -> prop_combo_box);


	QWidget * mapper_widget = new QWidget(this);
	QHBoxLayout * mapper_widget_layout = new QHBoxLayout(mapper_widget);
	mapper_widget_layout -> addWidget(new QLabel("Mapper",this));
	mapper_widget_layout -> addWidget(this -> mapper_combo_box);
	
	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);

	main_layout -> addWidget(shape_widget);
	main_layout -> addWidget(mapper_widget);
	main_layout -> addWidget(button_box);

	main_layout -> addStretch(1);

	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	
	this -> init();

}

void SelectMapperWindow::init(){
	
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		if (prop_combo_box -> count() == 1){
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("None"));
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Slopes"));
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Inertial Potential"));
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Body-fixed Potential"));

			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Inertial acceleration"));
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Body-fixed acceleration"));
			this -> mapper_combo_box -> insertItem(this -> mapper_combo_box -> count(),QString::fromStdString("Slope standard deviations"));

		}
	}

}


void SelectMapperWindow::accept(){

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();


	if (this -> mapper_combo_box -> count() > 0){

		int current_mapper_index = this -> mapper_combo_box -> currentIndex();

		if(current_mapper_index == 0){
			QDialog::accept();
		}

		vtkSmartPointer<vtkPolyData> shape = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_polydata();
		vtkSmartPointer<vtkPolyDataMapper> shape_mapper = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_mapper();
		vtkSmartPointer<vtkActor> shape_actor = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_actor();

		vtkSmartPointer<vtkFloatArray> data = nullptr;
		std::string unit;

		if (current_mapper_index == 0){
			shape_mapper -> ScalarVisibilityOff();
			this -> parent -> get_renderer() -> RemoveActor2D(wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_colorbar_actor());
			this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
			QDialog::accept();
		}
		else if (current_mapper_index == 1){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_slopes();
			unit = "(deg)";
		}
		else if (current_mapper_index == 2){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_inertial_potentials();
			unit = "(m^2/s^2)";

		}
		else if (current_mapper_index == 3){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_body_fixed_potentials();
			unit = "(m^2/s^2)";

		}
		else if (current_mapper_index == 4){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_inertial_acc_magnitudes();
			unit = "(m/s^2)";

		}
		else if (current_mapper_index == 5){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_body_fixed_acc_magnitudes();
			unit = "(m/s^2)";
			
		}
		else if (current_mapper_index == 6){
			data = wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_slope_sds();
			unit = "(deg)";
			
		}

		

		if (data == nullptr){
			QDialog::accept();
		}
		else{


			double valuesRange[2];
			data -> GetRange(valuesRange);

			if (std::isnan(valuesRange[0]) || std::isnan(valuesRange[1])){
				valuesRange[0] = 0;
				valuesRange[1] = 90;
			}

			double min_data = valuesRange[0];
			double max_data = valuesRange[1];

			shape -> GetCellData() -> SetScalars(data);
			shape -> Modified();

		  	// Create a lookup table to map cell data to colors
			vtkSmartPointer<vtkLookupTable> lut = vtkSmartPointer<vtkLookupTable>::New();
			lut -> SetNumberOfTableValues(shape -> GetNumberOfCells());
			lut -> SetRange(min_data,max_data);
			lut -> Build();

			shape_mapper -> SetScalarRange(min_data,max_data);
			shape_mapper -> SetLookupTable(lut);
			shape_mapper -> ScalarVisibilityOn();
			shape_mapper -> SetScalarModeToUseCellData();
			shape_mapper -> Update();

			shape_actor -> Modified();

			vtkSmartPointer<vtkScalarBarActor> scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
			scalarBar -> SetLookupTable(shape_mapper->GetLookupTable());
			scalarBar -> SetTitle((this -> mapper_combo_box -> currentText().toStdString()
				+ "\n"
				+ this -> prop_combo_box -> currentText().toStdString()
				+ "\n"
				+ unit).c_str()
			);

			scalarBar -> SetNumberOfLabels(4);


			// All of the other mappers/colorbars are removed
			for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
				if (it -> first != this -> prop_combo_box -> currentText().toStdString()){

					it -> second -> get_mapper() -> ScalarVisibilityOff();
					this -> parent -> get_renderer() -> RemoveActor2D(it -> second -> get_colorbar_actor());

				}
			}

			this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
			this -> parent -> get_renderer() -> RemoveActor2D(wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> get_colorbar_actor());
			this -> parent -> get_renderer() -> AddActor2D(scalarBar);
			wrapped_shape_data[this -> prop_combo_box -> currentText().toStdString()] -> set_colorbar_actor(scalarBar);
			this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

		}

	}

	QDialog::accept();

}


