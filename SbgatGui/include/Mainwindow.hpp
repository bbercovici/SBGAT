
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

#include <vtkSmartPointer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkAreaPicker.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>


#include "QVTKWidget.h"

#include <ShapeModelImporter.hpp>
#include <ShapeModel.hpp>

#include <map>



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
	Creates and display a vtkPolyData corresponding to the provided shape model
	@param shape_model Pointer to instantiated shape model
	*/
	void create_vtkpolydata_from_shape_model(SBGAT_CORE::ShapeModel * shape_model);



	QMenu * fileMenu;
	QMenu * ShapeModelMenu;
	QMenu * ViewMenu;

	vtkSmartPointer<vtkRenderer> renderer;
	vtkSmartPointer<vtkOrientationMarkerWidget> orientation_widget;
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor;

	std::shared_ptr<SBGAT_CORE::FrameGraph> frame_graph;

	std::map<std::string, std::shared_ptr<SBGAT_CORE::ShapeModel> > shape_models;


};

#endif