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

#include "ObsWindow.hpp"
#include <QMessageBox>
#include <QScrollArea>

#include <vtkSphereSource.h>
#include <vtkTriangleFilter.h>
#include "ShapePropertiesWidget.hpp"



using namespace SBGAT_GUI;

ObsWindow::ObsWindow(Mainwindow * parent) {

	this -> parent = parent;


	QGroupBox * target_group = new QGroupBox(tr("Shapes"));
		
	QScrollArea * scroll_area = new QScrollArea(this);
	// scroll_area -> setVerticalScrollBarPolicy( Qt::ScrollBarAlwaysOn );
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

