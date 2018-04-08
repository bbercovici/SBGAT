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


#ifndef HEADER_YORPWINDOW
#define HEADER_YORPWINDOW

#include <QMainWindow>
#include <QGroupBox>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>

#include <QLabel>
#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>

#include "Mainwindow.hpp"


namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class YORPWindow
\author Benjamin Bercovici
\date March, 2018
\brief YORPWindow class defining a window where a user can specificy the inputs to 
a YORP coefficients computation

\details This class inherits from QDialog and enables the
user to set a number of options used in the YORP coefficients calculation. The calculation 
can only take place if a spacecraft or a small body shape model is available and if an output directory 
has been specified
*/

	class YORPWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
		YORPWindow(Mainwindow * parent) ;


		private slots:

	/**
	Applies the currently seleted options to SbgatGUI
	*/
		void accept();

		/**
		Opens a dialog letting the user choose the output directory for the YORP coefficients
		*/

		void open_output_file_dialog();




	protected:

		void init();

		Mainwindow * parent;

		QComboBox * prop_combo_box;

		QComboBox * fourier_combo_box;
		QComboBox * bounces_combo_box;
		QComboBox * refine_combo_box;
		QComboBox * voxel_combo_box;

		QDialogButtonBox * button_box;

		QDoubleSpinBox * rho_sbox;
		QDoubleSpinBox * spec_sbox;

		QDoubleSpinBox * lambda_del_sbox;
		QDoubleSpinBox * delta_del_sbox;


		QPushButton * open_output_file_dialog_button;
		std::string output_path;

	};
}
#endif