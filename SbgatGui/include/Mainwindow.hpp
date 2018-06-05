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

/*!
@file Mainwindow.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief Stores definition of the Mainwindow class.
*/

#ifndef HEADER_MAINWINDOW
#define HEADER_MAINWINDOW

#include <QMainWindow>
#include <QDockWidget>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QWidget>
#include <QStatusBar>
#include <QPlainTextEdit>
#include <QTableWidget>

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPolyData.h>
#include <QVTKOpenGLWidget.h>


#include <map>
#include <chrono>
#include <sstream>
#include <memory>

#include "ModelDataWrapper.hpp"
#include <SBGATFrameGraph.hpp>


namespace SBGAT_GUI {

// Forward declaration of InteractorStyle
	class InteractorStyle;


	typedef std::map<std::string , std::shared_ptr<ModelDataWrapper> > DataMap;

/*!
@class Mainwindow
\author Benjamin Bercovici
\date July 22, 2017
\brief Mainwindow class. This is the main class of the SbgatGUI application

\details {Main class of the GUI as it hosts the VTK pipeline visualizer and
the actions/menus allowing the user to interact with the program data.
This class inherits from QMainWindow and hosts the QVTKOpenGLWidget where
all rendering and display tasks occur. Also exposes SbgatCore classes to the user
through the user interface layer brought by Qt.}
*/

	class Mainwindow : public QMainWindow {
		Q_OBJECT

	public:
	/**
	QVTKOpenGLWidget hosting the rendering pipeline.
	*/
		QVTKOpenGLWidget * qvtkWidget;


	/**
	Lateral dockwidget hosting a log console and
	shape selection options.
	*/
		QDockWidget * left_dockwidget;


	/**
	Lateral dockwidget hosting a log console and
	shape selection options.
	*/
		QDockWidget * right_dockwidget;

	/**
	Info bar providing SbgatGUI status and information
	about currently selected shape model.
	*/
		QStatusBar * status_bar;

	/**
	Read-only log console. Its content can be saved to a file.
	The console itself can be cleared from the menu bar.
	*/
		QPlainTextEdit * log_console;


	/**
	Shape selection table widget. Enables shape selection/ visibility toggle / erasing
	*/
		QTableWidget * prop_table;


	/**
	Setups the GUI and creates an instance of QVTK Widget.
	*/
		Mainwindow();

	/**
	Returns a pointer to the renderer associated with the window's QVTK widget.
	@return pointer to the vtkRenderer associated with the window's QVTK widget.
	*/
		vtkSmartPointer<vtkRenderer> get_renderer();

	/**
	Getter to wrapped shape data
	@return copy of wrapper shape data
	*/
		DataMap get_wrapped_shape_data() const;


	/**
	Returns to a pair storing the directory and vtkActor pointer to the current skybox
	@return a pair storing the directory and vtkActor pointer to the current skybox
	*/
		std::pair<std::string,vtkSmartPointer<vtkActor> > get_skybox_pair() const;

	/**
	Setter to the skybox's vtkActor
	@param skybox_actor Skybox actor
	*/
		void set_skybox_actor(vtkSmartPointer<vtkActor> skybox_actor);

	/**
	Setter to the skybox's directory
	@param skybox_dir Skybox directory
	*/
		void set_skybox_directory(std::string skybox_dir);



	/**
	Enable/Disables an action in the GUI
	@param enabled Status the targeted action will be set to
	@param action Pointer to action to enable/disable
	*/
		void set_action_status(bool enabled, QAction * action);


	// Actions


	/**
	When triggered, starts shape model saving action sequence.
	*/
		QAction * save_shape_action;

	/**
	When triggered, starts shape model loading action sequence.
	*/
		QAction * add_shape_action;

	

	/**
	When triggered, opens settings window.
	*/
		QAction * load_settings_window_action;

	/**
	When triggered, clears the log console.
	*/
		QAction * clear_console_action;

	/**
	When triggered, opens path dialog where a savepath is queried.
	If a valid path is provided, the content of the console is saved to this path
	*/
		QAction * save_console_action;

		/**
		When triggered, calls slot aligning the loaded small body shape
		with its bodycenter/principal axes
		*/
		QAction * align_shape_action;


	/**
	When triggered, show/hides the left lateral widget
	*/
		QAction * show_left_dockwidget_action;

	/**
	When triggered, show/hides the right lateral widget
	*/
		QAction * show_right_dockwidget_action;

	/**
	When triggered, opens a window enabling the user
	to align the displayed actors to a number of hook points
	corresponding to the center of mass of any of loaded props 
	*/
		QAction * open_alignment_window_action;



	/**
	When triggered, prints geometry measures of the selected prop to the console
	*/
		QAction * compute_geometric_measures_action;


	/**
	When triggered, opens settings window.
	*/
		QAction * open_settings_window_action;


	/**
	When triggered, opens a window letting the user compute the 
	Fourier coefficients of the YORP force/torque acting on a small body
	*/
		QAction * open_compute_yorp_window_action;


	/**
	When triggered, opens a window letting the user compute the spherical
	harmonics coeffients of the gravity acceleration
	*/
		QAction * open_compute_sharm_window_action;


		/**
	When triggered, opens camera properties window.
	*/
		QAction * open_rendering_properties_window_action;


		/**
		When triggered, opens the dialog window allowing
		one to compute simulated radar images over a targeted shape model
		*/

		QAction * open_radar_window_action;


	



		

		
		signals:

		/**
		Sends the signal that a prop has been added
		*/
		void prop_added_signal();

		/**
		Sends the signal that a prop has been removed
		*/
		void prop_removed_signal();



	// Slots
		private slots:

		




	/**
	Opens a widget allowing the user edit the camera properties
	*/
		void open_rendering_properties_window();




		/**
		Opens the dialog window allowing
		one to compute simulated radar images over a targeted shape model
		*/

		void open_radar_window();


	/**
	Removes the selected shape model from Sbgat by querying the name of the selected shape.
	Using this name to remove the shape, mappers, polydatas and actors associated with it.
	*/
		void remove_prop();


		/**
		Open a dialog window where the user can specify the parameters for the 
		computation of the Fourier expansion coefficients of the YORP force/torque 
		acting on a body
		*/
		void open_compute_yorp_window();


		/**
		Open a dialog window where the user can specify the parameters for the 
		computation of the spherical harmonics coefficients of the exterior gravity field
		*/
		void open_compute_sharm_window();

		/**
		Aligns the selected small body shape
		with its bodycenter/principal axes
		*/
		void align_shape();

		/**
		Computes and displays a number of geometry measures associated with the selected prop 
		*/
		void compute_geometric_measures();

		/**
	Shows/hides the selected prop from the lateral widget.
	@param row row index of calling cell
	@param col col index of calling cell (should be 1 for the slot to proceed, otherwise the call is ignored)
	*/
		void toggle_prop_visibility(int row, int col) ;


	/**
	Updates availability of GUI actions given latest model state.
	Ensures that all available actions
	are consistent with current Sbgat state.
	*/
		void update_actions_availability() ;


	/**
	Updates the GUI when a new shape model is selected from the lateral widget.
	*/
		void update_GUI_changed_prop();

	/**
	Open settings window.
	*/
		void open_settings_window();




	private:

		/**
		Add light of prescribed type to the renderer
		@param light_type determines the light type (0: scene light , 1: Head light , 2: camera light)
		*/
		// void add_light(int light_type);

		/**
		Initializes the rendering window and its props
		*/
		void init_rendering_window();

		/**	
		Initializes the right dockwidget. This widget holds information on loaded 
		props and the log console
		*/
		void init_right_dockwidget();


	/**
	Creates the GUI actions enabling the user to interact with the software, and connects them to the
	corresponding slots.
	*/
		void createActions();



	/**
	Adds a row in the table widget where shape models are listed,
	and fills it up to represent the newly loaded shape model on a new row.
	@param name name of new shape.
	*/
		void add_prop_to_table_widget(std::string name);


	/**
	Creates and populates the menu bar.
	*/
		void createMenus();

	/**
	Creates the GUI elements and places them in the main window.
	*/
		void setupUi();


	/**
	Load small body shape model stored in a .obj file. The shape model is stored in an instance of the ShapeModel
	class for subsequent operations. A vtkPolydata is also constructed for visualization purposes
	*/
		void add_shape();

		/**
	Save small body shape model to a .obj file.
	*/
		void save_shape();

	

	/**
	Clears the console.
	*/
		void clear_console() ;

	/**
	Saves the console content to a file.
	*/
		void save_console() ;

		QMenu * SettingsMenu;
		QMenu * SmallBodyMenu;
		QMenu * MeasuresMenu;
		QMenu * ObservationsMenu;
		QMenu * AnalysesMenu;
		QMenu * ConsoleMenu;
		QMenu * ResultsMenu;

		vtkSmartPointer<vtkRenderer> renderer;
		vtkSmartPointer<vtkOrientationMarkerWidget> orientation_widget;
		std::pair<std::string ,vtkSmartPointer<vtkActor> > skybox_pair;

		std::shared_ptr<SBGATFrameGraph> frame_graph;


		DataMap wrapped_shape_data;










	};
}
#endif