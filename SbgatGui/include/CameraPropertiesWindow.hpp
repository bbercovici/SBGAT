
#ifndef HEADER_CameraPropertiesWindow
#define HEADER_CameraPropertiesWindow

#include <QMainWindow>
#include <QVBoxLayout>
#include <QGroupBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QComboBox>
#include <vtkCamera.h>

#include <QVTKOpenGLWidget.h>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class MoveAlongTrajectoryWindow
\author Benjamin Bercovici
\date January 28, 2017
\brief CameraPropertiesWindow

\details 
*/

	class CameraPropertiesWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		CameraPropertiesWindow(Mainwindow * parent) ;


		public slots:
		/**
		Triggered after a prop has been removed from the main window. Ensures that a 
		deleted shape model/spacecraft can no longer be used 
		*/
		void prop_removed_slot();

		/**
		Triggered after a prop has been added to the main window. Will
		update the widget to reflect new prop
		*/
		void prop_added_slot();

		private slots:

		void change_focus();



		

	protected:

		Mainwindow * parent;
		void init();
		QComboBox * prop_combo_box;




	};
}
#endif