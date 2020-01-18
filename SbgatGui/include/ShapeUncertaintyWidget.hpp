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





#ifndef HEADER_SHAPE_UNCERTAINTY_WIDGET
#define HEADER_SHAPE_UNCERTAINTY_WIDGET

#include <QWidget>
#include <QDoubleSpinBox>
#include <QLabel>
#include <QTabWidget>
#include <QTableWidget>
#include "ObsWindow.hpp"

namespace SBGAT_GUI {

/*!
@class ShapeUncertaintyWidget
\author Benjamin Bercovici
\date April, 2019
\brief ShapeUncertaintyWidget class defining a widget where a user can specify the uncertainty
affecting the control points coordinates and rotation period of a considered shape

\details TODO
*/	

	class AnalysesWindow;
	

	class ShapeUncertaintyWidget : public QGroupBox {
		Q_OBJECT

	public:

	/**
	Creates the widget
	@param parent pointer to parent window.
	@param is_primary true if corresponding shape is a primary, false otherwise. Will add the option
	to set body density if true.
	@param title widget title
	*/
		ShapeUncertaintyWidget(AnalysesWindow * parent,std::string title) ;

		/**
		Returns user-set standard deviation in rotation periond
		@return standard deviation in rotation period (seconds)
		*/
		double get_period_sd() const;

		/**
		Returns density of corresponding body. If the body is not declared as 
		primary, this method will return 0
		@return density (kg/m^3)
		*/
		double get_density() const;


		/**
		Get index of current tab index
		@return active tab index
		*/
		int get_active_tab_index() const{return this -> tab_widget -> currentIndex();}


		/**
		Returns the string storing the input file from which the shape covariance is loaded
		@return covariance input file
		*/
		std::string get_covariance_input_file() const{return this -> covariance_input_file;}



		/**
		Returns the value of the global noise standard deviation (m)
		@return global noise standard deviation (m)
		*/
		double get_global_sigma() const {return this -> global_sigma_spinbox -> value();}
		
		/**
		Returns the value of the global noise correlation distance (m)
		@return global noise correlation distance (m)
		*/
		double get_global_correlation_distance() const {return this -> global_correlation_distance_spinbox -> value();}


		/**
		Get the user-specified number of covariance regularizations in the case of 
		global uncertainties
		@return number of times the shape covariance is regularized
		*/
		int get_global_covariance_regularization_number() const;
		/**
		Get the user-specified number of covariance regularizations in the case of 
		local uncertainties
		@return number of times the shape covariance is regularized
		*/
		int get_local_covariance_regularization_number() const;

		/**
		Returns state of save_global_covariance_to_file_checkbox checkbox
		*/
		bool get_save_global_covariance_to_file_checkbox() const{return this -> save_global_covariance_to_file_checkbox -> isChecked();}

		/**
		Returns state of save_local_covariance_to_file_checkbox checkbox
		*/
		bool get_save_local_covariance_to_file_checkbox() const{return this -> save_local_covariance_to_file_checkbox -> isChecked();}



		QTableWidget * local_uncertainty_table_widget;

		public slots:
		void clear();

		private slots:

		void select_covariance_input_file();
		void add_shape_uncertainty_region();
		void remove_shape_uncertainty_region();


	protected:

		AnalysesWindow * parent;


		QButtonGroup * uncertainty_type_button_group;
		QTabWidget * tab_widget;

		std::string covariance_input_file;
		QLabel * covariance_input_file_label;

		QDoubleSpinBox * global_sigma_spinbox;
		QDoubleSpinBox * global_correlation_distance_spinbox;

		QSpinBox * global_covariance_regularization_spin_box;
		QSpinBox * local_covariance_regularization_spin_box;

		QCheckBox *	save_global_covariance_to_file_checkbox;
		QCheckBox *	save_local_covariance_to_file_checkbox;



		QPushButton * covariance_input_file_button;
		QPushButton * add_region_button;
		QPushButton * remove_region_button;




		




	};
}
#endif