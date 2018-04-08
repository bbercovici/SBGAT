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

#include "RenderingPropertiesWindow.hpp"

#include <QVBoxLayout>
#include <QGroupBox>
#include <QPushButton>
#include <QCheckBox>

#include <QDialogButtonBox>
#include <QTableWidget>
#include <QLabel>
#include <vtkCamera.h>


#include <vtkLightCollection.h>
#include <vtkLight.h>
#include <vtkShadowMapPass.h>

#include <vtkSequencePass.h>
#include <vtkRenderPassCollection.h>
#include <vtkCameraPass.h>
#include <vtkOpenGLRenderer.h>
#include <vtkShadowMapBakerPass.h>


using namespace SBGAT_GUI;

RenderingPropertiesWindow::RenderingPropertiesWindow(Mainwindow * parent) : QDialog(parent,Qt::WindowStaysOnTopHint) {

	this -> parent = parent;
	this -> setWindowTitle("Rendering properties");

	QVBoxLayout * main_layout = new QVBoxLayout(this);
	QGroupBox * focus_prop_group = new QGroupBox(tr("Camera focus"));

	this -> prop_combo_box = new QComboBox (this);

	QLabel * focus_prop_label = new QLabel("Select prop to focus on",this);


	QGridLayout * focus_prop_layout = new QGridLayout(focus_prop_group);


	focus_prop_layout -> addWidget(focus_prop_label, 0, 0, 1, 1);
	focus_prop_layout -> addWidget(prop_combo_box, 0, 1, 1, 1);


	main_layout -> addWidget(focus_prop_group);
	

	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);


	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> prop_combo_box,SIGNAL(currentIndexChanged(int)),this,SLOT(change_focus()));
	
	this -> init();

}

void RenderingPropertiesWindow::init(){
	

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();
	
	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		
	}

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		
	}



}


void RenderingPropertiesWindow::prop_removed_slot(){

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

	for (int i = 0; i < this -> prop_combo_box -> count(); ++i){
		QString prop_name = this -> prop_combo_box -> itemText(i);
		if (wrapped_spacecraft_data.find(prop_name.toStdString()) == wrapped_spacecraft_data.end()&& wrapped_shape_data.find(prop_name.toStdString()) == wrapped_shape_data.end()){
			this -> prop_combo_box -> removeItem(i);
			break;
		}
	}
}







void RenderingPropertiesWindow::prop_added_slot(){
	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		if (this -> prop_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		if (this -> prop_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}



}

void RenderingPropertiesWindow::change_focus(){

	if (this -> prop_combo_box -> count() > 0 ){

		auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
		auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

		std::string current_prop_name = this -> prop_combo_box -> currentText().toStdString();

		vtkSmartPointer<vtkActor> prop_to_focus_on;

		if (wrapped_spacecraft_data.find(current_prop_name) != wrapped_spacecraft_data.end()){
			prop_to_focus_on = wrapped_spacecraft_data[current_prop_name] -> get_actor();
		}
		else{
			prop_to_focus_on = wrapped_shape_data[current_prop_name] -> get_actor();
		}



		vtkSmartPointer<vtkCamera> camera = this -> parent -> get_renderer() -> GetActiveCamera();
		camera -> SetFocalPoint(prop_to_focus_on -> GetPosition());

		this -> parent -> get_renderer() -> Modified();
		this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
	}


}

