
#ifndef HEADER_SELECTEDPOINTWIDGET
#define HEADER_SELECTEDPOINTWIDGET

#include <QDialog>
#include <QString>
#include <QTableWidget>
#include <QTableWidgetItem>
#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QLabel>
#include <QDialogButtonBox>
#include <QComboBox>
#include <QPushButton>
#include <QSlider>
#include <QLineEdit>


#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataNormals.h>

#include "interactor.h"
#include "mainwindow.h"

#include <set>
#include "vtkObjectFactory.h"


// Forward declaration of InteractorStyle
class InteractorStyle;

// Forward declaration of MainWindow
class MainWindow;


/**
Enum defining the state "transform direction" drop down list
*/
enum class TransformDirection {RADIAL, NORMAL};


/**
Enum defining the state "interpolation type" drop down list
*/
enum class InterpolationType {UNIFORM, LINEAR, PARABOLIC};


/**
Enum defining the state "transform selection" drop down list
*/
enum class TransformSelection {SELECTED, NCLOSEST};




/**
Declaration of the SelectedPointWidget class. SelectedPointWidget refers to the
widget displayed on the screen when the user selects at least one vertex of the
displayed shape model by means of the rectangular box selector. The widget that is then
displayed lists the IDs of the selected vertices, as well as a choice of possible geometric
transforms to be applied to them. A horizontal slider enalbes the user to choose
the magnitude of the transform
*/
class SelectedPointWidget : public QDialog {
	Q_OBJECT

public:
	/**
	Constructor.
	@param parent Pointer to parent widget (here, pointer to instance of MainWindow )
	*/

	SelectedPointWidget(QWidget *parent);

	/**
	Data Setter The pointer passed as arguments allow the widget to have access to the
	point properties
	@param interactor_style Pointer to the interactor accessing the shape model currently displayed
	*/
	void set_data(vtkSmartPointer<InteractorStyle> interactor_style);

	/**
	Highlights the cells corresponding to the selected points
	*/
	void highlight_selected_cells();

	/**
	Computes the normals at the selected points
	*/
	void compute_selected_points_normals();



	/**
	Reset the widget by clearing the memory used by the different vtk container members

	*/
	void reset();

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

	// Slots
private slots:
	/**
	Displays the current slider position
	@param pos Slider position
	*/
	void show_new_slider_pos(int pos);

	/**
	Redraw the shape model based on the slider position
	@param pos Slider position
	*/
	void update_view(int pos);

	/**
	Sets the position slider value via the QLineEdit field
	*/
	void set_new_slider_pos();

	void accept();
	void reject();

	/**
	Sets the states of transform_direction_list
	*/
	void set_transform_direction(const int);


	/**
	Sets the states of interpolation_type_list
	*/
	void set_interpolation_type(const int);

	/**
	Sets the states of transform_selection_list
	*/
	void set_transform_selection(const int);


private:
	void createActions();
	void createMenus();
	bool * widget_is_open;

	/**
	Loops over the list of all actors created by this and remove each of those from the rendering window
	*/
	void remove_selected_points_actor();
	/**
	Populates the QTableWidget table with the relevant data
	*/
	void populate_vertex_table();

	/**
	Constructs the vtkPolyData corresponding to 1) the selected cells 2) the rest of the cells in the
	currently displayed shape model. This allows the rendering window to represent both sets of cells
	separately
	*/
	void compute_cell_blobs();

	/**
	Computes the memory space used by all the class members derived from vtkAbstractArray
	@return size Int representing the memory used by class members derived from
	vtkAbstractArray in MbiBytes
	*/
	float get_actual_memory_size();


	QStringList labels;

	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	vtkSmartPointer<vtkPolyData> all_points_polydata;
	vtkSmartPointer<vtkPolyData> selected_cells_polydata;
	vtkSmartPointer<vtkPolyData> unselected_cells_polydata;



	vtkSmartPointer<vtkPoints> selected_points;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	vtkSmartPointer<vtkCellArray> selected_polys_cell_array;
	vtkSmartPointer<vtkCellArray> unselected_polys_cell_array;
	vtkSmartPointer<vtkIdTypeArray> selected_polys_ids;
	vtkSmartPointer<vtkIdTypeArray> unselected_polys_ids;
	vtkSmartPointer<vtkPoints> new_selected_points_coordinates;

	vtkSmartPointer<vtkIdTypeArray> polys_ids;
	vtkSmartPointer<vtkIdTypeArray> visible_points_global_ids_from_local_index;

	vtkSmartPointer<vtkIdList> cell_ids;
	vtkDataArray * selected_points_normals;

	vtkSmartPointer<vtkPolyDataNormals> selected_polydata_normals_filter;

	TransformDirection transform_direction;
	InterpolationType interpolation_type;
	TransformSelection transform_selection;


	MainWindow * mainwindow;
};




#endif

