
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
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkAreaPicker.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkArrowSource.h>

#include "QVTKWidget.h"
#include "osxHelper.h"


// forward declaration of InteractorStyle
class InteractorStyle;

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
	Returns a pointer to the window interactor
	*/
	vtkSmartPointer<vtkRenderWindowInteractor> get_render_window_interactor();



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

	/**
	Shows/hides facet normals
	*/

	QAction * show_facet_normals_action;


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
	std::vector<vtkSmartPointer<vtkActor> > normal_actors;



};

#endif