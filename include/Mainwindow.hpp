
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

#include "QVTKWidget.h"
#include "osxHelper.h"


#include "SelectedPointWidget.hpp"
#include "Interactor.hpp"

// forward declaration of InteractorStyle
class InteractorStyle;
class SelectedPointWidget;

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
	Sets the visibility of all the actors owned by this
	@param visibility Boolean setting the visibility of the actors owned by this
	*/
	void set_actors_visibility(bool visibility);

	/**
		Enable/Disables an action in the GUI
		@param enabled Status the targeted action will be set to
		@param action Pointer to action to enable/disable
		*/
	void set_action_status(bool enabled, QAction * action);

	QAction * openAct;
	QAction * saveAct;
	QAction * selectPointAct;
	QAction * shapeColorAct;
	QAction * resetAct;
	QAction * backgroundColorAct;
	QAction * vertexVisibilityAct;


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
	Enables/Disables vertex visibility on the displayed shape model
	*/
	void change_vertex_visibility();


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

	QMenu * fileMenu;
	QMenu * OperationMenu;
	QMenu * ViewMenu;
	
	vtkSmartPointer<vtkRenderer> renderer;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;



};

#endif