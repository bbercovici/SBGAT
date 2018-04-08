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