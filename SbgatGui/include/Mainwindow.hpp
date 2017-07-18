
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
#include <QCoreApplication>

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


#include "QVTKWidget.h"

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>

#include <map>
#include <chrono>
#include <sstream>



// forward declaration of InteractorStyle
class InteractorStyle;

/**
Declaration of the Mainwindow Class. Main class of the GUI as it hosts the VTK pipeline visualizer and
the actions/menus allowing the user to interact with the program data.
*/

class Mainwindow : public QMainWindow {
	Q_OBJECT

public:
	// QVTKWidget
	QVTKWidget * qvtkWidget;

	// Docks
	QDockWidget * lateral_dockwidget;

	// Status Bar
	QStatusBar * status_bar;

	// Log console
	QPlainTextEdit * log_console;



	/**
	Returns a pointer to the vtkRenderer associated with the window's QVTK widget
	@return Pointer to the vtkRenderer associated with the window's QVTK widget
	*/
	vtkSmartPointer<vtkRenderer> get_renderer();

	/**
	Constructor. Setups the GUI and creates an instance of QVTK Widget
	*/
	Mainwindow();


	/**
	Closes any opened lateral dockwidget
	*/
	void close_lateral_dockwidget();

	/**
	Returns a pointer to the window interactor
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> get_render_window_interactor();


	/**
	Enable/Disables an action in the GUI
	@param enabled Status the targeted action will be set to
	@param action Pointer to action to enable/disable
	*/
	void set_action_status(bool enabled, QAction * action);





	// Actions

	/**
	Open load model
	*/
	QAction * load_shape_model_action;

	/**
	Save shape model
	*/

	QAction * save_action;

	/**
	Modify shape model
	*/
	QAction * modify_shape_action;

	/**
	Change shape color
	*/
	QAction * set_shape_color_action;

	/**
	Clear currently open shape model
	*/
	QAction * clear_all_action;

	/**
	Set background color
	*/
	QAction * set_background_color_action;

	/**
	Open ComputePGMWidget
	*/
	QAction * open_ComputePGMWidget_action;

	/**
	Open ShapeInfoWidget
	*/
	QAction * open_ShapeInfoWidget_action;

	/**
	Shows/hides facet normals
	*/

	QAction * show_facet_normals_action;

	/**
	Shows/hides lateral dockwidget
	*/
	QAction * show_lateral_dockwidget_action;


	/**
	Clears the log console
	*/
	QAction * clear_console_action;

	/**
	Saves the log console to a file
	*/
	QAction * save_console_action;

	/**
	Saves the log console to a file
	*/
	QAction * print_inertia_action;

	/**
	Evaluates the polyhedron gravity model at the specified point in the
	principal body frame
	*/
	QAction * compute_pgm_acceleration_action;


	/**
	Evaluates the polyhedron gravity model at center of each facet, evaluated in the
	principal body frame
	*/
	QAction * compute_global_pgm_acceleration_action;

	/**
	Evaluates gravitational slope at the center of each facet
	*/
	QAction * compute_gravity_slopes_action;

	/**
	Evaluates gravitational slope at the center of each facet
	*/
	QAction * toggle_grav_slopes_visibility_action;


	// Slots
private slots:

	/**
	Sets the background color to that chosen by the user
	*/
	void set_background_color();



private:
	/**
	Creates the GUI actions enabling the user to interact with the software, and connects them to the
	appropriate slots
	*/
	void createActions();

	/**
	Creates and populates the menu bar
	*/
	void createMenus();

	/**
	Creates the GUI elements and places them in the main window
	*/
	void setupUi();


	/**
	Load shape model stored in a .obj file. The shape model is stored in an instance of the ShapeModel
	class for subsequent operations. A vtkPolydata is also constructed for visualization purposes
	*/
	void load_shape_model();


	/**
	Enables/disables visibility of gravity slopes, provided that they have been computed
	*/
	void toggle_grav_slopes_visibility();

	/**
	Transfer results of gravity slopes computation into vtk
	*/
	void update_vtk_slopes() ;

	/**
	Prints dimensionless principal inertia tensor of the target (assuming constant density)
	to the log console
	*/
	void print_inertia() const ;


	/**
	Computes the polyhedron gravity model acceleration at the specified point in the
	shape's principal body frame
	*/
	void compute_pgm_acceleration();


	/**
	Computes the polyhedron gravity model acceleration at the center of each facet,
	expressed in the shape's principal body frame
	*/
	void compute_global_pgm_acceleration() ;


	/**
	Computes the gravity slopes at the center of each facet.
	WARNING: the gravity acceleration at each facet must have been completed first
	*/
	void compute_gravity_slopes() ;



	/**
	Creates and display a vtkPolyData corresponding to the provided shape model
	@param shape_model Pointer to instantiated shape model
	@param name name under which the polydata should be saved (same as the shape model)
	*/
	void create_vtkpolydata_from_shape_model(SBGAT_CORE::ShapeModel * shape_model,
	        std::string name);

	/**
	Shows/hides lateral dockwidget
	*/
	void show_lateral_dockwidget();

	/**
	Clears the console
	*/
	void clear_console() ;

	/**
	Saves the console content to a file
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
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

	std::shared_ptr<SBGAT_CORE::FrameGraph> frame_graph;

	std::map<std::string, std::shared_ptr<SBGAT_CORE::ShapeModel> > shape_models;
	std::map<std::string, vtkSmartPointer<vtkPolyData> > polydatas;
	std::map<std::string, vtkSmartPointer<vtkPolyDataMapper> > mappers;
	std::map<std::string, vtkSmartPointer<vtkActor> > actors;




};

#endif