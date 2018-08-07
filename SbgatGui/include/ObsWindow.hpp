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
\brief ObsWindow class defining a window where a user can generate simulated observations involving (optionally)binary asteroids.


\details The position and attitude states of the asteroids can either be generated under keplerian dynamics (resp. fixed-spin) regime, or provided 
from an external input file. 
Note:
- Should the position of a secondary be generated under the keplerian dynamics assumption, the center of mass of the primary will 
remain at (0,0,0) with no velocity (in violation of the actual barycentric motion of the two bodies). 
- The attitude is parametrized in terms of a set of Modified Rodrigues Parameters sigma = (sigma_1,sigma_2,sigma_3)
- Position and attitude files should be header-less, whitespace-separated, time-sorted like so :

t_0 a_0 b_0 c_0 d_0 e_0 f_0
t_1 a_1 b_1 c_1 d_1 e_1 f_1
...........................
t_N a_N b_N c_N d_N e_N f_N

where t_i is the state epoch expressed in seconds, (a_i, b_i, c_i) the position cartesian coordinates (resp. unitless MRP)  expressed in meters ( resp. unitless), 
(d_i, e_i, f_i) the velocity (resp. angular velocity) expressed in meters/second (resp. expressed in rad/s). 



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
		Loads time history of 6x1 state into the provided containers. 
		Expects to read file line-by-line with 7 values on each line
		(time, followed by the 6 state components) with each value separated by 
		a whitespace
		@param filepath path to file
		@param state_vec container to hold each vector from the state history
		@param time_vec container to hold each time in the state history
		*/
		static void load_state_from_file(const std::string & filepath,
			std::vector<arma::vec> & state_vec,
			std::vector<double> & time_vec);
		
		void get_inputs_from_GUI(std::vector<double> & imaging_times,
			std::vector< std::vector<arma::vec> > & positions_vec, 
			std::vector< std::vector<arma::vec> > & velocities_vec,
			std::vector< std::vector<arma::vec> > & mrps_vec,
			std::vector< std::vector<arma::vec> > & omegas_vec);



		void add_state_history(ShapePropertiesWidget * shape_properties_widget,
			const std::vector<double> & imaging_times,
			std::vector< std::vector<arma::vec> > & positions_vec, 
			std::vector< std::vector<arma::vec> > & velocities_vec,
			std::vector< std::vector<arma::vec> > & mrps_vec,
			std::vector< std::vector<arma::vec> > & omegas_vec) const;

		virtual void init();

		/**
		Piecewise linear interpolation of provided state history
		@param time interpolation time
		@param state_from_file vector of states , each of them being 6x1
		@param state_time_from_file vector of state epochs
		@param is_attitude_state true if the interpolated state is an attitude state (sigma_1,sigma_2,sigma_3,omega_1,omega_2,omega_3).
		Will ensure that the interpolator behaves properly around MRP shadow set switching, if any
		@return interpolated state
		*/
		static arma::vec interpolate(const double & time,
			const std::vector<arma::vec> & state_from_file,
			const std::vector<double> & state_time_from_file,
			bool is_attitude_state);

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