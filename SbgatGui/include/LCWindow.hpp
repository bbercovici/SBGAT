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


#ifndef HEADER_LC_WINDOW
#define HEADER_LC_WINDOW

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
#include <SBGATObsLightcurve.hpp>



namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class LCWindow
\author Benjamin Bercovici
\date March, 2018
\brief LCWindow class defining a window where a user can generate emulated radar 
data simulating the output of a range/range-rate doppler radar

\details TODO
*/

	class LCWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
		LCWindow(Mainwindow * parent) ;


		private slots:


		/**
		Collect lightcurve with specified inputs
		*/	
		void collect_observations();


		/**
		Save lightcurve data to png and raw .txt file
		*/

		void save_observations();

		/**
		Opens the visualizer to view lightcurve
		*/
		void open_visualizer();


		/**
		Update displayed value of phase angle
		*/
		void update_phase_angle();



	protected:

		void init();

		Mainwindow * parent;

		QComboBox * prop_combo_box;

		QDialogButtonBox * button_box;

		QDoubleSpinBox * spin_raan_sbox;
		QDoubleSpinBox * spin_inc_sbox;

		QDoubleSpinBox * observer_az_sbox;
		QDoubleSpinBox * observer_el_sbox;

		QDoubleSpinBox * sun_az_sbox;
		QDoubleSpinBox * sun_el_sbox;

		QDoubleSpinBox * rotation_period_sbox;
		QDoubleSpinBox * imaging_period_sbox;

		QSpinBox * N_samples_sbox;
		QSpinBox * N_images_sbox;

		QCheckBox * penalize_incidence_box;


		QPushButton * save_observations_button;
		QPushButton * collect_observations_button;
		QPushButton * open_visualizer_button;

		QLabel * phase_angle_label;
		QLineEdit * phase_angle_qldt;

		std::string output_path;

		vtkSmartPointer<SBGATObsLightcurve> lc;

		std::vector<std::array<double, 2> > measurements;

	};
}
#endif