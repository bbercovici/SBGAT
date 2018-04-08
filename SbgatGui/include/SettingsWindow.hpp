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