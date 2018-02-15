
#ifndef HEADER_SETTINGSWINDOW
#define HEADER_SETTINGSWINDOW

#include <QMainWindow>
#include <QGroupBox>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class SettingsWindow
\author Benjamin Bercovici
\date July 22, 2017
\brief SettingsWindow class enabling preferences setting.

\details This class inherits from QDialog and enables the
user to set a number of options. Default options values will eventually be
saved to a hidden configuration file loaded anytime SbgatGui starts up.
*/

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

	/**
	Opens up QFileDialog where the user can select a directory
	storing 6 .png images forming a skybox cubic texture 
	*/
		void open_skybox_directory_dialog();

	/**	
	Displays the selected background type options
	*/
		void show_selected_background_type(int index);



	protected:

		Mainwindow * parent;

		QComboBox * aa_frames_combo_box;
		QComboBox * background_type_combo_box;

		QPushButton * open_color_dialog_button;
		QPushButton * open_skybox_directory_dialog_button;

		QCheckBox * use_gradient_checkbox;
		QGroupBox * background_color_group;
		QGroupBox * skybox_group ;



		double rgb_button[3];


	};
}
#endif