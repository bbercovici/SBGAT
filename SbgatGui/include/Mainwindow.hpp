/*!
@file Mainwindow.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief Stores definition of the Mainwindow class.
*/

#ifndef HEADER_MAINWINDOW
#define HEADER_MAINWINDOW

#include <QMainWindow>
#include <QFileDialog>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QDockWidget>
#include <QtWidgets/QAction>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QColorDialog>
#include <QColor>
#include <QtWidgets/QWidget>
#include <QStatusBar>
#include <QInputDialog>
#include <QPlainTextEdit>
#include <QTextStream>
#include <QMessageBox>
#include <QRegularExpression>
#include <QStringList>
#include <QTableWidget>
#include <QPushButton>
#include <QHeaderView>
#include <QThread>

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkAreaPicker.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkScalarBarActor.h>
#include <vtkActor2DCollection.h>
#include <vtkLookupTable.h>
#include <vtkTextProperty.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkParametricSpline.h>
#include <vtkParametricFunctionSource.h>
#include <vtkOBJReader.h>
#include <vtkCenterOfMass.h>
#include <QVTKOpenGLWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>

#include <map>
#include <chrono>
#include <sstream>

#include "SettingsWindow.hpp"
#include "ModelDataWrapper.hpp"
#include "MoveAlongTrajectoryWindow.hpp"
#include "CameraPropertiesWindow.hpp"


#include "Worker.hpp"



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
		QDockWidget * lateral_dockwidget;

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
	Getter to wrapped trajectory data
	@return copy of wrapper trajectory data
	*/
	DataMap get_wrapped_trajectory_data() const;

	/**
	Getter to wrapped attitude data
	@return copy of attitude shape data
	*/
	DataMap get_wrapped_attitude_data() const;

	/**
	Getter to wrapped spacecraft data
	@return copy of wrapper spacecraft data
	*/
	DataMap get_wrapped_spacecraft_data() const;

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
	When triggered, starts shape model loading action sequence.
	*/
		QAction * load_small_body_action;

	/**
	When triggered, starts trajectory loading action sequence.
	*/
		QAction * load_trajectory_action;


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
	When triggered, show/hides the lateral widget
	*/
		QAction * show_lateral_dockwidget_action;

	/**
	When triggered, prints geometry measures of the selected prop to the console
	*/
		QAction * compute_geometry_measures_action;


	/**
	When triggered, starts polyhedron gravity model computation sequence by
	querying point coordinates in the
	principal body frame and bulk density of attracting body.
	*/
		QAction * compute_pgm_acceleration_action;

	/**
	When triggered, opens settings window.
	*/
		QAction * open_settings_window_action;


		/**
	When triggered, opens camera properties window.
	*/
		QAction * open_camera_properties_window_action;

	/**
	When triggered, starts global polyhedron gravity model accelerations evaluation. Queries the
	bulk density of attracting body before computing pgm accelerations at the center of each facet.
	*/
		QAction * compute_global_pgm_acceleration_action;


	/**
	When triggered, starts global polyhedron gravity model potential evaluation. Queries the
	bulk density of attracting body before computing pgm potential at the center of each facet.
	*/
		QAction * compute_global_pgm_potential_action;

	/**
	When triggered, starts evaluation sequence of gravitational slope at the center of each facet.
	Will query orientation of spin axis along with spin rate. Only available if the
	global PGM accelerations of the selected shape have already been computed
	*/
		QAction * compute_grav_slopes_action;


	/**
	When triggered, opens up a widget where visibility of
	computed gravitational slopes can be turned on/off.
	*/
		QAction * show_grav_slopes_action;


	/**
	When triggered, opens up a widget where visibility of
	computed gravitational potentials can be turned on/off.
	*/
		QAction * show_global_pgm_pot_action;




	/**
	When triggered, opens the dialog window
	allowing one to load a spacecraft shape model
	*/
		QAction * load_spacecraft_action;

	/**
	When triggered, opens the dialog window
	allowing one to move a previously loaded spacecraft along a previously loaded trajectory
	*/
		QAction * move_along_traj_action;
		
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
	Opens a widget allowing the user to associate a trajectory to a spacecraft 
	and move the spacecraft along it
	*/
		void open_move_along_traj_window();


	/**
	Opens a widget allowing the user edit the camera properties
	*/
		void open_camera_properties_window();


	/**
	Removes the selected shape model from Sbgat by querying the name of the selected shape.
	Using this name to remove the shape, mappers, polydatas and actors associated with it.
	*/
		void remove_prop();


		/**
		Computes and displays a number of geometry measures associated with the selected prop 
		*/
		void compute_geometry_measures();

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


	/**
	Transfer results of gravity slopes computation into VTK.
	*/
		void update_vtk_slopes() ;

	/**
	Transfer results of gravity potential computation into VTK.
	*/
		void update_vtk_potentials();


	private:

	/**
	Creates the GUI actions enabling the user to interact with the software, and connects them to the
	corresponding slots.
	*/
		void createActions();


	/**
	Removes props associated with display of results
	@param name name of shape model associated with the removal of the visual props
	@param remove_all true if all props should be removed. Otherwise only those corresponding to the
	name in argument will be removed.
	*/
		void remove_results_visual_props(std::string name, bool remove_all) ;


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
		void load_small_body();

	/**
	Load spacecraft shape model stored in a .obj file. A vtkPolydata is also constructed for visualization purposes
	*/
		void load_spacecraft();






	/**
	Load x/y/z trajectory expressed in a small body's body-fixed frame
	*/
		void load_trajectory();


	/**
	Show gravitational slopes by choosing one shape model for which gravitational slopes are available.
	*/
		void show_grav_slopes();


	/**
	Show gravitational potentials by choosing one shape model for which gravitational potentials are available.
	*/
		void show_global_pgm_pot();



	/**
	Prints dimensionless principal inertia tensor of the active shape (assuming constant density)
	to the log console.
	*/
		void compute_inertia() ;


	/**
	Prints volume of active shape (m^3) to the log console.
	*/
		void compute_volume() ;


	/**
	Prints center of mass coordinates shape (m) to the log console.
	*/
		void compute_center_of_mass() ;


	/**
	Prints surface of active shape (m^2) to the log console.
	*/
		void compute_surface_area() ;


	/**
	Computes the polyhedron gravity model acceleration at the specified point in the
	shape's principal body frame.
	*/
		void compute_pgm_acceleration();




	/**
	Computes the polyhedron gravity model acceleration at the center of each facet,
	expressed in the shape's principal body frame.
	*/
		void compute_global_pgm_acceleration() ;


	/**
	Computes the polyhedron gravity model potential at the center of each facet
	*/
		void compute_global_pgm_potential() ;


	/**
	Computes the gravity slopes at the center of each facet.
	The gravity acceleration at each facet must have been computed first, so this
	method is not accessible before it is the case for the currently selected
	shape model.
	*/
		void compute_gravity_slopes() ;





	/**
	Creates and display a vtkPolyData corresponding to the provided shape model.
	Also stores the associated vtkPolyDataMapper and vtkActor
	@param model_data pointer to the ModelDataWrapper housing the data
	related to the shape model being created
	*/
		void create_vtkpolydata_from_shape_model(std::shared_ptr<ModelDataWrapper> model_data);

	/**
	Shows/hides lateral dockwidget.
	*/
		void show_lateral_dockwidget();

	/**
	Clears the console.
	*/
		void clear_console() ;

	/**
	Saves the console content to a file.
	*/
		void save_console() ;

		QMenu * SmallBodyMenu;
		QMenu * SpacecraftMenu;
		QMenu * MeasuresMenu;
		QMenu * TrajectoryMenu;
		QMenu * ViewMenu;
		QMenu * DynamicAnalysesMenu;
		QMenu * ConsoleMenu;
		QMenu * ResultsMenu;

		vtkSmartPointer<vtkRenderer> renderer;
		vtkSmartPointer<vtkOrientationMarkerWidget> orientation_widget;
		std::pair<std::string ,vtkSmartPointer<vtkActor> > skybox_pair;

		std::shared_ptr<SBGAT_CORE::FrameGraph> frame_graph;


		DataMap wrapped_shape_data;
		DataMap wrapped_trajectory_data;
		DataMap wrapped_attitude_data;
		DataMap wrapped_spacecraft_data;









	};
}
#endif