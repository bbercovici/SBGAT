
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
user to set a number of options. Default options values will eventually be
saved to a hidden configuration file loaded anytime SbgatGui starts up.
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