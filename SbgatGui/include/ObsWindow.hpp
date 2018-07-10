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


#ifndef HEADER_OBS_WINDOW
#define HEADER_OBS_WINDOW

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
#include <SBGATObs.hpp>



namespace SBGAT_GUI {

	class Mainwindow;
	class ShapePropertiesWidget;

/*!
@class ObsWindow
\author Benjamin Bercovici
\date March, 2018
\brief ObsWindow class defining a window where a user can generate emulated radar 
data simulating the output of a range/range-rate doppler radar

\details TODO
*/

	class ObsWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
		ObsWindow(Mainwindow * parent) ;

		private slots:

		void changed_secondary_box(int index);


	protected:

		/**


		*/
		void get_inputs_from_GUI(std::vector<double> & imaging_times,
			std::vector< std::vector<arma::vec> > & positions_vec, 
			std::vector< std::vector<arma::vec> > & velocities_vec,
			std::vector< std::vector<arma::vec> > & mrps_vec,
			std::vector< std::vector<arma::vec> > & omegas_vec);


		virtual void init();

		Mainwindow * parent;

		QComboBox * primary_prop_combo_box;
		QComboBox * secondary_prop_combo_box;


		QDialogButtonBox * button_box;
		QDoubleSpinBox * imaging_period_sbox;


		QSpinBox * N_samples_sbox;
		QSpinBox * N_images_sbox;

		QCheckBox * penalize_incidence_box;

		QGroupBox * obs_specific_group;
		QPushButton * save_observations_button;
		QPushButton * open_visualizer_button;
		QPushButton * collect_observations_button;
		ShapePropertiesWidget * primary_shape_properties_widget;
		ShapePropertiesWidget * secondary_shape_properties_widget;

		vtkSmartPointer<SBGATObs> observation_filter;

	};
}
#endif