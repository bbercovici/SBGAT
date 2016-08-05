
#ifndef HEADER_SELECTEDPOINTWIDGET
#define HEADER_SELECTEDPOINTWIDGET

#include <QDialog>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QDialogButtonBox>
#include <QComboBox>
#include <QPushButton>
#include <QSlider>

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkIdTypeArray.h>

#include "interactor.h"
#include "mainwindow.h"

#include <set>

// Forward declaration of InteractorStyle
class InteractorStyle;

// Forward declaration of MainWindow
class MainWindow;


/**
Declaration of the SelectedPointWidget class. SelectedPointWidget refers to the
widget displayed on the screen when the user selects at least one vertex of the
displayed shape model by means of the rectangular box selector. The widget that is then
displayed lists the IDs of the selected vertices, as well as a choice of possible geometric
transforms to be applied to them
NOTE: for now, the only possible transform (homothetic transform) is hardcoded
*/
class SelectedPointWidget : public QDialog {
	Q_OBJECT

public:
	/**
	Constructor. 
	*/
	
	SelectedPointWidget();

	/**
	Data Setter The pointer passed as arguments allow the widget to have access to the
	point properties
	@param points_polydata Pointer to the vtkPolyData storing the vertices of the shape model
	@param selected_points_polydata Pointer to the vtkPolyData storing the selected vertices
	of the displayed shape model
	*/
	void set_data(vtkSmartPointer<InteractorStyle> interactor_style);

	/**
	Highlights the cells corresponding to the selected points
	*/
	void highlight_selected_cells();

	QTableWidget * table;
	QHBoxLayout * slider_layout;
	QVBoxLayout * main_layout;
	QSlider * slider;
	QLineEdit * slider_value;
	QDialogButtonBox * button_box;
	QWidget * slider_holder_widget;
	QLabel * transform_direction_title;
	QLabel * interpolation_type_title;
	QLabel * transform_selection_title;
	QLabel * slider_title;


	QComboBox * transform_direction_list;
	QComboBox * interpolation_type_list;
	QComboBox * transform_selection_list;


	QPushButton * button_show_vertex_table;



	// Slots
private slots:
	void show_new_slider_pos(int pos);
	void update_view(int pos);

	void set_new_slider_pos();
	void show_vertex_table();
	void accept();
	void reject();


private:
	void createActions();
	void createMenus();
	bool * widget_is_open;
	void remove_selected_points_actor();
	/**
	Populates the QTableWidget table with the relevant data
	*/
	void populate_vertex_table();

	/**
	Returns a pointer to the structured polydata representing the selected facets/vertices
	*/
	vtkSmartPointer<vtkPolyData> get_selected_blob_polydata();



	QStringList labels;
	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	vtkSmartPointer<vtkPolyData> points_polydata;
	vtkSmartPointer<vtkPolyData> selected_cells_polydata;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	MainWindow * mainwindow;
};




#endif

