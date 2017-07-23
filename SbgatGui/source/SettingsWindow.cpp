#include "SettingsWindow.hpp"

using namespace SBGAT_GUI;

SettingsWindow::SettingsWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Preferences");

	QVBoxLayout * settings_layout = new QVBoxLayout(this);
	QGroupBox * antialiasing_group = new QGroupBox(tr("Antialiasing"));
	QGroupBox * background_color_group = new QGroupBox(tr("Background color"));
	QGridLayout * antialiasing_group_layout = new QGridLayout(antialiasing_group);
	QGridLayout * background_color_group_layout = new QGridLayout(background_color_group);


	this ->  open_color_dialog_button = new QPushButton(this);
	this -> use_gradient_checkbox = new QCheckBox(this);

	QLabel * label_samples = new QLabel("AA frames", this);
	QLabel * label_background_color = new QLabel("Choose color", this);
	QLabel * label_use_gradient = new QLabel("Use gradient background", this);


	this -> aa_frames_combo_box = new QComboBox (this);

	for (unsigned int i = 0; i < 13 ; ++i) {
		aa_frames_combo_box -> insertItem(i, QString::number(i));
	}

	// Getting the current antialiasing settings
	aa_frames_combo_box -> setCurrentIndex(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetAAFrames());
	antialiasing_group_layout -> addWidget(label_samples, 0, 0, 1, 1);
	antialiasing_group_layout -> addWidget(aa_frames_combo_box, 0, 1, 1, 1);


	// Creating the button box
	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
	        | QDialogButtonBox::Cancel);


	// Background color
	background_color_group_layout -> addWidget(label_background_color, 0, 0, 1, 1);
	background_color_group_layout -> addWidget(open_color_dialog_button, 0, 1, 1, 1);

	background_color_group_layout -> addWidget(label_use_gradient, 1, 0, 1, 1);
	background_color_group_layout -> addWidget(use_gradient_checkbox, 1, 1, 1, 1);

	this -> use_gradient_checkbox -> setChecked(this -> parent -> get_renderer() -> GetGradientBackground ());


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

	// Adding elements to the main layout
	settings_layout -> addWidget(antialiasing_group);
	settings_layout -> addWidget(background_color_group);
	settings_layout -> addWidget(button_box);

	// Connecting buttons signals to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> open_color_dialog_button, SIGNAL(clicked(bool)), this, SLOT(open_color_dialog()));




}


void SettingsWindow::accept() {

	// Set new AA settings
	this -> parent -> qvtkWidget -> GetRenderWindow() -> SetAAFrames(int(this -> aa_frames_combo_box -> currentText().toDouble()));

	// Updating the parent window with the new color settings and rendering
	this -> parent -> get_renderer() -> SetBackground(
	    this -> rgb_button[0] / 255.,
	    this -> rgb_button[1] / 255.,
	    this -> rgb_button[2] / 255.);
	this -> parent -> get_renderer() -> SetGradientBackground(this -> use_gradient_checkbox -> isChecked());
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


	QDialog::accept();

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
