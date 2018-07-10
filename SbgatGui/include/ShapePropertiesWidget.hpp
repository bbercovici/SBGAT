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





#ifndef HEADER_SHAPE_PROPERTIES_WIDGET
#define HEADER_SHAPE_PROPERTIES_WIDGET

#include <QWidget>
#include <QDoubleSpinBox>
#include <QLabel>

#include "ObsWindow.hpp"


namespace SBGAT_GUI {

/*!
@class ShapePropertiesWidget
\author Benjamin Bercovici
\date March, 2018
\brief ShapePropertiesWidget class defining a widget where a user can specify the values taken 
by the physical parameters of a shape

\details TODO
*/

	class ShapePropertiesWidget : public QGroupBox {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	@param is_primary true if corresponding shape is a primary, false otherwise. Will add the option
	to set body density if true.
	@param title widget title
	*/
		ShapePropertiesWidget(ObsWindow * parent,bool is_primary,std::string title) ;

		/**
		Returns direction of spin vector in inertial frame
		@return spin vector (unit vector)
		*/
		arma::vec get_spin() const;

		/**
		Returns rotation periond
		@return rotation period (seconds)
		*/
		double get_period() const;

		/**
		Returns density of corresponding body. If the body is not declared as 
		primary, this method will return 0
		@return density (kg/m^3)
		*/
		double get_density() const;


		/**
		Returns a vector of keplerian orbital elements
		@return vector of orbital elements sma (m) ,e ,i (rad) ,Omega (rad) ,omega (rad) ,M0 (rad)
		*/
		arma::vec get_orbital_elements() const;


		

		/**
		Exclusive push button which, if pressed, indicates that the position of the body is keplerian
		*/
		QRadioButton * position_from_keplerian_button ;

		/**
		Exclusive push button which, if pressed, indicates that the position of the body is to be loaded from a file
		*/
		QRadioButton * position_from_file_button ;

		/**
		Exclusive push button which, if pressed, indicates that the attitude of the body derives from a fixed-spin regime
		*/
		QRadioButton * attitude_from_simple_spin_button ;


		/**
		Exclusive push button which, if pressed, indicates that the attitude of the body is to be loaded from a file
		*/
		QRadioButton * attitude_from_file_button ;

		private slots:

		void toggle_position_visibility();
		void toggle_attitude_visibility();



	protected:

		void init();

		QLabel * spin_raan_label;
		QLabel * spin_inc_label;
		QLabel * period_label;



		QDoubleSpinBox * spin_raan_sbox;
		QDoubleSpinBox * spin_inc_sbox;
		QDoubleSpinBox * period_sbox;



		QDoubleSpinBox * sma_sbox;
		QDoubleSpinBox * ecc_sbox;
		QDoubleSpinBox * inc_sbox;
		QDoubleSpinBox * Omega_sbox;
		QDoubleSpinBox * omega_sbox;
		QDoubleSpinBox * M0_sbox;



		QDoubleSpinBox * density_sbox = nullptr;
		QLabel * density_label = nullptr;

		

		QWidget * position_keplerian_widget;
		QWidget * attitude_simple_spin_widget;


		QWidget * position_from_file_widget;
		QWidget * attitude_from_file_widget;

		QGroupBox * position_box;
		QGroupBox * attitude_box;

		

	};
}
#endif