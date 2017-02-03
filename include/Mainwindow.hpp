
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
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QWidget>
#include <QStatusBar>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointSource.h>
#include <vtkCleanPolyData.h>
#include <vtkDelaunay3D.h>
#include <vtkCellArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkOBJReader.h>
#include <vtkOBJExporter.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkCamera.h>
#include <vtkLine.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPolyDataNormals.h>

#include "QVTKWidget.h"
#include "osxHelper.h"


#include "ModifyAreaWidget.hpp"
#include "InteractorStyle.hpp"
#include "ComputePGMWidget.hpp"
#include "ShapeInfoWidget.hpp"
#include "SetInputScalingWidget.hpp"
#include "Asteroid.hpp"

#include "vtkOBJWriter.h"

// forward declaration of InteractorStyle
class InteractorStyle;
class ModifyAreaWidget;

/**
Declaration of the Mainwindow Class. Main class of the GUI as it hosts the VTK pipeline visualizer and
the actions/menus allowing the user to interact with the program data.
*/

class Mainwindow : public QMainWindow {
	Q_OBJECT

public:
	// Widgets
	QVTKWidget * qvtkWidget;

	// Docks
	QDockWidget * lateral_dockwidget;

	// Status Bar
	QStatusBar * status_bar;

	/**
	Returns a pointer to the vtkRenderer associated with the window's QVTK widget
	@return Pointer to the vtkRenderer associated with the window's QVTK widget
	*/
	vtkSmartPointer<vtkRenderer> get_renderer();

	/** Constructor. Setups the GUI and creates an instance of QVTK Widget
	*/
	Mainwindow();

	/**
	Closes any opened lateral dockwidget
	*/
	void close_lateral_dockwidget();

	/**
	Return a copy of the vector storing the pointers to the vtk actors currently displayed
	@return Vector of vtkActor *
	*/
	std::vector<vtkSmartPointer<vtkActor> > get_actor_vector();


	/**
	Returns a pointer to the window interactor
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> get_render_window_interactor();

	/**
	Sets the visibility of all the actors owned by this
	@param visibility Boolean setting the visibility of the actors owned by this
	*/
	void set_actors_visibility(bool visibility);

	/**
	Returns a pointer to the asteroid currently loaded
	@return Pointer to the asteroid 
	*/
	Asteroid * get_asteroid();

	/**
	Sets the scaling factor applied to the 
	input data to have it expressed in meters
	*/
	void set_scaling_factor(double scaling_factor);

	/**
		Enable/Disables an action in the GUI
		@param enabled Status the targeted action will be set to
		@param action Pointer to action to enable/disable
		*/
	void set_action_status(bool enabled, QAction * action);

	/**
	Open load model 
	*/
	QAction * load_action;

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


	// Slots
private slots:
	/**
	Opens a shape model
	*/
	void open();

	/**
	Saves the vertex data of the shape model currently displayed to an $prefix.obj file with $prefix chosen by the user
	This function also generates $prefix.mtl file which contains material/lightning information
	*/
	void save();

	/**
	Allows the interactor to grab props by
	setting its style mode to INTERACTOR_IS_SELECT
	It is rigourosly equivalent to pressing the "r" key
	*/
	void select();

	/**
	Sets the shape color to that chosen by the user. A shape model must be displayed so that
	the corresponding action is enabled
	*/
	void set_shape_color();

	/**
	Sets the background color to that chosen by the user
	*/
	void set_background_color();

	/**
	Opens a widget displaying information on the 
	currently opened shape model
	*/
	void open_shape_info_widget();

	
	/**
	Opens a widget allowing the user to compute the polyhedron gravity model of the displayed shape model
	*/
	void open_compute_pgm_widget();


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
	Takes in a pointer to a vtkPolyData representing a shape model,
	generates the convex hull corresponding to this dataset and
	displays both on the main window.
	*/
	void load_obj(vtkSmartPointer<vtkPolyData> read_polydata_without_id);

	/**
	Creates the GUI elements and places them in the main window
	*/
	void setupUi();

	/**
	Remove all the actors owned by this instance of mainwindow
	*/
	void remove_actors();

	/**
	Clear all windows and loaded shape models
	*/
	void clear_all();

	/**
	Scaling factor applied to the 
	input data to have it expressed in meters
	*/
	double scaling_factor;

	QMenu * fileMenu;
	QMenu * ShapeModelMenu;
	QMenu * ViewMenu;

	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkOrientationMarkerWidget> widget;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

	Asteroid * asteroid;



};

#endif