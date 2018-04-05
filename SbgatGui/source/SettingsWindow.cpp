/** MIT License

Copyright (c) 2018 Benjamin Bercovici

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

#include "SettingsWindow.hpp"


#include <QLabel>
#include <QVBoxLayout>
#include <QColorDialog>
#include <QColor>
#include <QDialogButtonBox>

#include <vtkTexture.h>
#include <vtkSkybox.h>
#include <vtkImageFlip.h>
#include <vtkPNGReader.h>
#include <QVTKOpenGLWidget.h>


using namespace SBGAT_GUI;

SettingsWindow::SettingsWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Preferences");

	QVBoxLayout * settings_layout = new QVBoxLayout(this);
	QGroupBox * antialiasing_group = new QGroupBox(tr("Antialiasing"));
	QGroupBox * background_type_group = new QGroupBox(tr("Background type"));
	this ->  background_color_group = new QGroupBox(tr("Background color"));
	this ->  skybox_group = new QGroupBox(tr("Skybox directory"));

	QGridLayout * antialiasing_group_layout = new QGridLayout(antialiasing_group);
	QGridLayout * background_type_group_layout = new QGridLayout(background_type_group);

	QGridLayout * background_color_group_layout = new QGridLayout(background_color_group);
	QGridLayout * skybox_group_layout = new QGridLayout(skybox_group);



	this ->  open_color_dialog_button = new QPushButton(this);
	this ->  open_skybox_directory_dialog_button = new QPushButton(this);

	this -> use_gradient_checkbox = new QCheckBox(this);

	QLabel * label_samples = new QLabel("AA frames", this);
	QLabel * label_background_color = new QLabel("Choose color", this);
	QLabel * label_background_type = new QLabel("Choose background type", this);
	QLabel * label_use_gradient = new QLabel("Use gradient background", this);
	QLabel * label_skybox_directory = new QLabel("Choose skybox directory", this);


	this -> aa_frames_combo_box = new QComboBox (this);
	this ->  background_type_combo_box = new QComboBox (this);
	for (unsigned int i = 0; i < 13 ; ++i) {
		this -> aa_frames_combo_box -> insertItem(i, QString::number(i));
	}
	this -> background_type_combo_box -> insertItem(0,"Color");
	this -> background_type_combo_box -> insertItem(1,"Skybox");

	// Getting the current antialiasing settings
	aa_frames_combo_box -> setCurrentIndex(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetAAFrames());
	antialiasing_group_layout -> addWidget(label_samples, 0, 0, 1, 1);
	antialiasing_group_layout -> addWidget(aa_frames_combo_box, 0, 1, 1, 1);

	// Background type group
	background_type_group_layout -> addWidget(label_background_type, 0, 0, 1, 1);
	background_type_group_layout -> addWidget(this -> background_type_combo_box, 0, 1, 1, 1);

	// Background color group
	background_color_group_layout -> addWidget(label_background_color, 0, 0, 1, 1);
	background_color_group_layout -> addWidget(open_color_dialog_button, 0, 1, 1, 1);
	background_color_group_layout -> addWidget(label_use_gradient, 1, 0, 1, 1);
	background_color_group_layout -> addWidget(use_gradient_checkbox, 1, 1, 1, 1);
	this -> use_gradient_checkbox -> setChecked(this -> parent -> get_renderer() -> GetGradientBackground ());

	// Skybox group
	skybox_group_layout -> addWidget(label_skybox_directory, 0, 0, 1, 1);
	skybox_group_layout -> addWidget(open_skybox_directory_dialog_button, 0, 1, 1, 1);

	// Getting the current color from the renderer in the parent window
	this -> parent -> get_renderer() -> GetBackground(this -> rgb_button);
	this -> rgb_button[0] = 255 * this -> rgb_button[0];
	this -> rgb_button[1] = 255 * this -> rgb_button[1];
	this -> rgb_button[2] = 255 * this -> rgb_button[2];
	QString color_string = QString::fromStdString("background-color: rgb(" +
		std::to_string(rgb_button[0])
		+  ","
		+ std::to_string(rgb_button[1])
		+ "," +
		std::to_string(rgb_button[2])
		+ ");" +
		"border-style: outset;" +
		"border-width: 2px;" +
		"min-width: 5em;" +
		"border-color: beige;" +
		"border-radius: 10px;");

	this -> open_color_dialog_button -> setStyleSheet(color_string);

	// Creating the button box
	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);

	// Adding elements to the main layout
	settings_layout -> addWidget(antialiasing_group);
	settings_layout -> addWidget(background_type_group);
	settings_layout -> addWidget(background_color_group);
	settings_layout -> addWidget(skybox_group);
	settings_layout -> addWidget(button_box);

	// Connecting buttons signals to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> open_color_dialog_button, SIGNAL(clicked(bool)), this, SLOT(open_color_dialog()));
	connect(this -> open_skybox_directory_dialog_button, SIGNAL(clicked(bool)), this, SLOT(open_skybox_directory_dialog()));
	connect(this -> background_type_combo_box,SIGNAL(currentIndexChanged(int)),this,SLOT(show_selected_background_type(int)));


	// Initializing 
	this -> background_color_group -> hide();
	this -> skybox_group -> hide();
	background_type_combo_box -> setCurrentIndex(-1);
	if (this -> parent -> get_skybox_pair(). first == ""){
		this -> background_type_combo_box -> setCurrentIndex(0);
	}
	else{
		this -> background_type_combo_box -> setCurrentIndex(1);
		this -> open_skybox_directory_dialog_button -> setText(QString::fromStdString(this -> parent -> get_skybox_pair(). first));
	}




}


void SettingsWindow::show_selected_background_type(int index){

	if (index == 0){
		this -> background_color_group -> show();
		this -> skybox_group -> hide();
	}
	else{
		this -> background_color_group -> hide();
		this -> skybox_group -> show();

	}

}

void SettingsWindow::accept() {

	// Set new AA settings
	this -> parent -> qvtkWidget -> GetRenderWindow() -> SetAAFrames(int(this -> aa_frames_combo_box -> currentText().toDouble()));

	if (open_skybox_directory_dialog_button -> text() != "" && this -> background_type_combo_box -> currentIndex() == 1
		&& open_skybox_directory_dialog_button -> text().toStdString() != this -> parent -> get_skybox_pair().first){
		std::string skybox_path = open_skybox_directory_dialog_button -> text().toStdString();
		const std::string fpath[] = {
			skybox_path + "/pX.png",
			skybox_path + "/nX.png",
			skybox_path + "/pY.png",
			skybox_path + "/nY.png",
			skybox_path + "/pZ.png",
			skybox_path + "/nZ.png",
		};

		vtkSmartPointer<vtkTexture> texture = vtkSmartPointer<vtkTexture>::New() ;
		texture->CubeMapOn();
		texture->InterpolateOn();
		texture->RepeatOff();
		texture->EdgeClampOn();

		texture->MipmapOn();

		for (int i = 0; i < 6; i++){
			vtkNew<vtkPNGReader> imgReader;
			imgReader->SetFileName(fpath[i].c_str());
			vtkNew<vtkImageFlip> flip;
			flip->SetInputConnection(imgReader->GetOutputPort());
			flip->SetFilteredAxis(1); 
			texture->SetInputConnection(i, flip->GetOutputPort(0));
		}


  // Create the skybox
		vtkSmartPointer<vtkSkybox> skybox_actor = vtkSmartPointer<vtkSkybox>::New();

		skybox_actor -> SetTexture(texture);

		this -> parent -> get_renderer() -> AddActor(skybox_actor);
		this -> parent -> set_skybox_directory(skybox_path);
		this -> parent -> set_skybox_actor(skybox_actor);
	}
	else if(this -> background_type_combo_box -> currentIndex() == 1 && 
		open_skybox_directory_dialog_button -> text().toStdString() == this -> parent -> get_skybox_pair().first){
		//same skybox, do nothing
	}
	else{
		if (this -> parent -> get_skybox_pair(). first != ""){
			this -> parent -> get_renderer() -> RemoveActor(this -> parent -> get_skybox_pair(). second);
			this -> parent -> set_skybox_directory("");
			this -> parent -> set_skybox_actor(nullptr);
		}


	}


	// Updating the parent window with the new color settings and rendering
	this -> parent -> get_renderer() -> SetBackground(
		this -> rgb_button[0] / 255.,
		this -> rgb_button[1] / 255.,
		this -> rgb_button[2] / 255.);
	this -> parent -> get_renderer() -> SetGradientBackground(this -> use_gradient_checkbox -> isChecked());
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


	QDialog::accept();

}

void SettingsWindow::open_skybox_directory_dialog(){

	QString dir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),
		"~",
		QFileDialog::ShowDirsOnly
		| QFileDialog::DontResolveSymlinks);



	QStringList cube_sides ;
	QStringList found_cube_sides ;
	cube_sides << "pX.png" << "nX.png" << "pY.png" << "nY.png"  << "pZ.png"  << "nZ.png";


	QDir currentDir(dir);
	foreach (const QString &match, currentDir.entryList(cube_sides, QDir::Files | QDir::NoSymLinks))
	found_cube_sides.append(match);

	if (found_cube_sides.size() == 6){
		open_skybox_directory_dialog_button -> setText(dir);
	}
	else{
		QMessageBox msgBox;
		msgBox.setText("This folder does not contain the expected cube textures. They are either missing or labelled incorrectly");
		msgBox.exec();
		open_skybox_directory_dialog_button -> setText("");


	}


}


void SettingsWindow::open_color_dialog() {

	QColor default_color(int(this -> rgb_button[0]),
		int(this -> rgb_button[1]),
		int(this -> rgb_button[2]));



	QColor color = QColorDialog::getColor(default_color);

	if (color.isValid()) {

		int red;
		int green;
		int blue;

		color.getRgb(&red, &green, &blue);
		QString color_string = QString::fromStdString("background-color: rgb(" +
			std::to_string(red) +  ","
			+ std::to_string(green) + "," +
			std::to_string(blue) + ");" +
			"border-style: outset;" +
			"border-width: 2px;" +
			"min-width: 5em;" +
			"border-color: beige;" +
			"border-radius: 10px;");

		this -> open_color_dialog_button -> setStyleSheet(color_string);
		this -> rgb_button[0] = red;
		this -> rgb_button[1] = green;
		this -> rgb_button[2] = blue;



	}





}
