/**
Mainwindow.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief Mainwindow class. This is the main class of the SbgatGUI application

\details This class inherits from QMainWindow and hosts the QVTKOpenGLWidget where
all rendering and display tasks occur. Also exposes SbgatCore classes to the user
through the user interface layer brought by Qt.
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

#include <QVTKOpenGLWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>

#include <map>
#include <chrono>
#include <sstream>

#include "SettingsWindow.hpp"
#include "ModelDataWrapper.hpp"


// Forward declaration of InteractorStyle
class InteractorStyle;


namespace SBGAT_GUI {


/**
Main class of the GUI as it hosts the VTK pipeline visualizer and
the actions/menus allowing the user to interact with the program data.
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
	QTableWidget * shape_table;


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
	Enable/Disables an action in the GUI
	@param enabled Status the targeted action will be set to
	@param action Pointer to action to enable/disable
	*/
	void set_action_status(bool enabled, QAction * action);


	// Actions

	/**
	When triggered, starts shape model loading action sequence.
	*/
	QAction * load_shape_model_action;


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
	When triggered, prints the shape inertia to the log console
	*/
	QAction * print_inertia_action;


	/**
	When triggered, prints the shape volume to the log console
	*/
	QAction * print_volume_action;

	/**
	When triggered, prints the shape surface area to the log console
	*/
	QAction * print_surface_action;

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
	When triggered, starts global polyhedron gravity model accelerations evaluation. Queries the
	bulk density of attracting body before computing pgm accelerations at the center of each facet.
	*/
	QAction * compute_global_pgm_acceleration_action;

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


	// Slots
private slots:


	/**
	Removes the selected shape model from Sbgat by querying the name of the selected shape.
	Using this name to remove the shape, mappers, polydatas and actors associated with it.
	*/
	void remove_shape();

	/**
	Shows/hides the selected shape model from the lateral widget.
	@param row row index of calling cell
	@param col col index of calling cell (should be 1 for the slot to proceed, otherwise the call is ignored)
	*/
	void toggle_shape_visibility(int row, int col) ;


	/**
	Updates availability of GUI actions given latest model state.
	Ensures that all available actions
	are consistent with current Sbgat state.
	*/
	void update_actions_availability() ;


	/**
	Updates the GUI when a new shape model is selected from the lateral widget.
	*/
	void update_GUI_changed_shape_model();

	/**
	Open settings window.
	*/
	void open_settings_window();


private:

	/**
	Creates the GUI actions enabling the user to interact with the software, and connects them to the
	corresponding slots.
	*/
	void createActions();


	/**
	Removes props associated with display of gravity slopes.
	@param name name of shape model associated with the removal of the visual props
	@param remove_all true if all props should be removed. Otherwise only those corresponding to the
	name in argument will be removed.
	*/
	void remove_grav_slopes_props(std::string name, bool remove_all) ;


	/**
	Adds a row in the table widget where shape models are listed,
	and fills it up to represent the newly loaded shape model on a new row.
	@param name name of new shape.
	*/
	void add_shape_to_table_widget(std::string name);


	/**
	Creates and populates the menu bar.
	*/
	void createMenus();

	/**
	Creates the GUI elements and places them in the main window.
	*/
	void setupUi();


	/**
	Load shape model stored in a .obj file. The shape model is stored in an instance of the ShapeModel
	class for subsequent operations. A vtkPolydata is also constructed for visualization purposes
	*/
	void load_shape_model();


	/**
	Show gravitational slopes by choosing one shape model for which gravitational slopes are available.
	*/
	void show_grav_slopes();

	/**
	Transfer results of gravity slopes computation into VTK.
	*/
	void update_vtk_slopes() ;

	/**
	Prints dimensionless principal inertia tensor of the active shape (assuming constant density)
	to the log console.
	*/
	void print_inertia() ;


	/**
	Prints volume of active shape (m^3) to the log console.
	*/
	void print_volume() ;


	/**
	Prints surface of active shape (m^2) to the log console.
	*/
	void print_surface() ;


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
	Computes the gravity slopes at the center of each facet.
	The gravity acceleration at each facet must have been completed first, so this
	method is not accessible before it is the case for the currently selected
	shape model.
	*/
	void compute_gravity_slopes() ;



	/**
	Creates and display a vtkPolyData corresponding to the provided shape model
	@param shape_model pointer to instantiated shape model
	@param name name under which the polydata should be saved
	*/
	void create_vtkpolydata_from_shape_model(SBGAT_CORE::ShapeModel * shape_model,
	        std::string name);

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

	QMenu * FileMenu;
	QMenu * ShapeMenu;
	QMenu * ViewMenu;
	QMenu * DynamicAnalysesMenu;
	QMenu * ConsoleMenu;
	QMenu * ResultsMenu;



	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkOrientationMarkerWidget> orientation_widget;

	std::shared_ptr<SBGAT_CORE::FrameGraph> frame_graph;

	std::map<std::string, std::shared_ptr<SBGAT_CORE::ShapeModel> > shape_models;
	std::map<std::string, vtkSmartPointer<vtkPolyData> > polydatas;
	std::map<std::string, vtkSmartPointer<vtkPolyDataMapper> > mappers;
	std::map<std::string, vtkSmartPointer<vtkActor> > actors;

	std::map<std::string, bool > consistent_global_accelerations;
	std::map<std::string, bool > consistent_grav_slopes;

	std::map<std::string , std::shared_ptr<ModelDataWrapper> > wrapped_data;







};
}
#endif