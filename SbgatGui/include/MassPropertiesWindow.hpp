/** MIT License

Copyright (c) 2019 Benjamin Bercovici and Jay McMahon

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


#ifndef HEADER_MASSPROPERTIES_WINDOW
#define HEADER_MASSPROPERTIES_WINDOW

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
#include "AnalysesWindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;
	class ShapePropertiesWidget;
	class ShapeUncertaintyWidget;

/*!
@class MassPropertiesWindow
\author Benjamin Bercovici
\date October, 2018
\brief MassPropertiesWindow class defining a window where a user 
evaluate shape models mass properties (volume, surface area, inertia,...)
and associated uncertainty
*/

	class MassPropertiesWindow : public AnalysesWindow {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		MassPropertiesWindow(Mainwindow * parent) ;
		
		private slots:
		
		void compute_mass_properties();

		

	protected:


		virtual void init();


		QDialogButtonBox * button_box;
		QCheckBox * save_to_file_checkbox;

		QPushButton * compute_mass_properties_button;
		ShapeUncertaintyWidget * primary_shape_uncertainty_widget;

		std::string output_path;
		
	};
}
#endif