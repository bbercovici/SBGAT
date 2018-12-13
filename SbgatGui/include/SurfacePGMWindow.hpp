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


#ifndef HEADER_SURFACEPGM_WINDOW
#define HEADER_SURFACEPGM_WINDOW

#include <QMainWindow>
#include <QGroupBox>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>

#include <QLabel>
#include <QLineEdit>

#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QDoubleSpinBox>
#include <QSpinBox>
#include <QRadioButton>


#include "Mainwindow.hpp"


namespace SBGAT_GUI {

	class Mainwindow;
	class ShapePropertiesWidget;

/*!
@class SurfacePGMWindow
\author Benjamin Bercovici
\date October, 2018
\brief SurfacePGMWindow class defining a window where a user 
evaluate the Polyhedron Gravity Model of a shape model
\details Enables computation of surface potential, inertial accelerations, 
body-frame accelerations and surface slopes
*/

	class SurfacePGMWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		SurfacePGMWindow(Mainwindow * parent) ;

		private slots:
		
		void compute_surface_pgm();
		void load_surface_pgm();

		/**
		Opens a dialog letting the user choose the output file for the computed surface polyhedron gravity model
		*/

		void open_output_file_dialog();




	protected:


		virtual void init();

		Mainwindow * parent;

		QComboBox * primary_prop_combo_box;

		QDialogButtonBox * button_box;
		QPushButton * open_output_file_dialog_button;

		QPushButton * compute_surface_pgm_button;
		QPushButton * load_surface_pgm_button;

		ShapePropertiesWidget * primary_shape_properties_widget;

		std::string output_path;
		
	};
}
#endif