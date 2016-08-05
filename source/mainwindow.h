
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
#include "QVTKWidget.h"


#include "newShapeModelDialog.h"
#include "selectedPointWidget.h"
#include "Interactor.h"

// forward declaration of InteractorStyle
class InteractorStyle;
class SelectedPointWidget;

/**
Declaration of the MainWindow Class. Main class of the GUI as it hosts the VTK pipeline visualizer and
the actions/menus allowing the user to interact with the program data.
*/

class MainWindow : public QMainWindow {
	Q_OBJECT

public:
	// Widgets
	QVTKWidget * qvtkWidget;
    SelectedPointWidget * pc_editing_widget;
    
	// Docks
	QDockWidget * selected_point_dockwidget;
	
	/**
	Returns a pointer to the vtkRenderer associated with the window's QVTK widget
	@return Pointer to the vtkRenderer associated with the window's QVTK widget
	*/
	vtkSmartPointer<vtkRenderer> get_renderer();

	/** Constructor. Setups the GUI and creates an instance of QVTK Widget
	*/
	MainWindow();

	/**
	Destructor
	*/
	~MainWindow();
    
    /**
     Pointer to boolean indicating whether the selection widget is already open
     */
    bool * selection_widget_is_open;
    
    

	// Slots
private slots:
	/**
	Opens a .vtp file storing vertex data using the path provided by the user
	*/
	void open();

	/**
	Saves the vertex data of the shape model currently displayed to a .vtp file whose full name is chosen by the user
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
	Generates a new shape model with N points normally distributed in space.
	*/
	void new_shape_model();

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
	Takes in a pointer to a vtkPolyData representing a point cloud,
	generates the convex hull corresponding to this dataset and
	displays both on the main window.
	*/
	void load_pc(vtkSmartPointer<vtkPolyData> read_polydata_without_id);

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
	Enable/Disables an action
	@param enabled Status the targeted action will be set to
	@param menu_name String storing the name of the menu containing the targeted action
	@param action_name String storing the name of the targeted action
	*/
	void set_action_status(bool enabled, const std::string & menu_name, const std::string & action_name);

	
	QMenu * fileMenu;
	QMenu * OperationMenu;
	QMenu * ViewMenu;
	QAction * openAct;
	QAction * saveAct;
	QAction * newShapeModelAct;
	QAction * selectPointAct;
	QAction * shapeColorAct;
	QAction * backgroundColorAct;
	QAction * vertexVisibilityAct;
    QColorDialog * palette;


	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<InteractorStyle> style;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;



};

#endif