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


#ifndef HEADER_RADAR_WINDOW
#define HEADER_RADAR_WINDOW

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
#include <QSpinBox>
#include <QRadioButton>

#include "Mainwindow.hpp"
#include <SBGATObsRadar.hpp>



namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class RadarWindow
\author Benjamin Bercovici
\date March, 2018
\brief RadarWindow class defining a window where a user can generate emulated radar 
data simulating the output of a range/range-rate doppler radar

\details TODO
*/

	class RadarWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
		RadarWindow(Mainwindow * parent) ;


		private slots:


		/**
		Collect radar observations with specified inputs
		*/	
		void collect_observations();


		/**
		Bins raw range/range-rate measurements to binned imaged
		*/

		void bin_observations();

		/**
		Save binned radar observations to PNG files
		*/

		void save_observations();


		/**
		Opens the visualizer to view computed radar images
		*/
		void open_visualizer();




	protected:

		void init();

		Mainwindow * parent;

		QComboBox * prop_combo_box;

		QDialogButtonBox * button_box;

		QDoubleSpinBox * r_bin_sbox;
		QDoubleSpinBox * rr_bin_sbox;
		QDoubleSpinBox * spin_raan_sbox;
		QDoubleSpinBox * spin_inc_sbox;

		QDoubleSpinBox * radar_az_sbox;
		QDoubleSpinBox * radar_el_sbox;

		QDoubleSpinBox * rotation_period_sbox;
		QDoubleSpinBox * imaging_period_sbox;

		QSpinBox * N_samples_sbox;
		QSpinBox * N_images_sbox;

		QPushButton * save_observations_button;
		QPushButton * collect_observations_button;
		QPushButton * bin_observations_button;
		QPushButton * open_visualizer_button;


		std::string output_path;

		SBGATMeasurementsSequence measurement_sequence;
		vtkSmartPointer<SBGATObsRadar> radar;

	};
}
#endif