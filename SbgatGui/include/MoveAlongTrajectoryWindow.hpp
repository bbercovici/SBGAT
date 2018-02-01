
#ifndef HEADER_MOVEALONGTRAJECTORYWINDOW
#define HEADER_MOVEALONGTRAJECTORYWINDOW

#include <QMainWindow>
#include <QGroupBox>
#include <QSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QDialog>


#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class MoveAlongTrajectoryWindow
\author Benjamin Bercovici
\date January 28, 2017
\brief MoveAlongTrajectoryWindow class enabling spacecraft motion along precomputed trajectories

\details 
*/

	class MoveAlongTrajectoryWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		MoveAlongTrajectoryWindow(Mainwindow * parent) ;


		public slots:
		/**
		Triggered after a prop has been removed from the main window. Ensures that a 
		deleted trajectory/spacecraft can no longer be used 
		*/
		void prop_removed_slot();

		/**
		Triggered after a prop has been added to the main window. Will
		update the widget to reflect new prop
		*/
		void prop_added_slot();

		private slots:
		/**
		Updates the position of the spacecraft after either the slider or
		the spin box have been modified. Makes the slider and the spin
		box consistent with each other
		*/
		void update_position();

		/**
		Updates the window elements when a new trajectory is selected
		*/
		void changed_trajectory();


		/**
		Updates the camera focus setting
		*/
		void toggle_camera_focus(int state);


	protected:

		Mainwindow * parent;
		void init();


		QComboBox * trajectory_combo_box;
		QComboBox * spacecraft_combo_box;
		QComboBox * camera_orientation_combo_box;
		QCheckBox * camera_focus_check;

		QGroupBox * position_group ;
		QSlider * pos_slider;
		QSpinBox * pos_spinbox;



	};
}
#endif