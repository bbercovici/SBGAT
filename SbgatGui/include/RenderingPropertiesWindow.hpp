
#ifndef HEADER_RenderingPropertiesWindow
#define HEADER_RenderingPropertiesWindow

#include <QMainWindow>
#include <QComboBox>
#include <QVTKOpenGLWidget.h>
#include <QDialog.h>
#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class RenderingPropertiesWindow
\author Benjamin Bercovici
\date January 28, 2017
\brief RenderingPropertiesWindow
\details A window where the user can set the camera focus, add/remove lights and enable/disable shadows
*/

	class RenderingPropertiesWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the window
	@param parent pointer to parent window.
	*/
		RenderingPropertiesWindow(Mainwindow * parent) ;


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
		// void add_light();
		// void remove_light();
		// void enable_mutual_shadows(int state);




	protected:

		Mainwindow * parent;
		void init();
		QComboBox * prop_combo_box;

		// QComboBox * current_light_combo_box;
		// QPushButton * remove_light_button;
		// QComboBox * new_light_combo_box ;

		// void make_light_box_consistent();
		





	};
}
#endif