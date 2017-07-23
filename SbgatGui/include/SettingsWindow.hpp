/**
SettingsWindow.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief SettingsWindow class enabling preferences setting.

\details This class inherits from QDialog and enables the
user to set a number of options. Default options values will eventually be
saved to a hidden configuration file loaded anytime SbgatGui starts up
*/


#ifndef HEADER_SETTINGSWINDOW
#define HEADER_SETTINGSWINDOW

#include <QMainWindow>
#include <QVBoxLayout>
#include <QColorDialog>
#include <QColor>
#include <QGroupBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QComboBox>
#include <QCheckBox>

#include <QVTKOpenGLWidget.h>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {

class Mainwindow;

class SettingsWindow : public QDialog {
	Q_OBJECT

public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
	SettingsWindow(Mainwindow * parent) ;


private slots:
	
	/**
	Applies the currently seleted options to SbgatGUI
	*/
	void accept();

	/**
	Opens up QColorPalette where a background color can be selected.
	*/
	void open_color_dialog();


protected:

	Mainwindow * parent;

	QComboBox * aa_frames_combo_box;
	QPushButton * open_color_dialog_button;
	QCheckBox * use_gradient_checkbox;


	double rgb_button[3];


};
}
#endif