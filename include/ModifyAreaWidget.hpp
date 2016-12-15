
#ifndef HEADER_ModifyAreaWidget
#define HEADER_ModifyAreaWidget

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
#include <vtkCellData.h>
#include <vtkIdTypeArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkKdTreePointLocator.h>
#include <vtkPoints.h>
#include <vtkCenterOfMass.h>

#include "vtkObjectFactory.h"

#include "InteractorStyle.hpp"
#include "Mainwindow.hpp"

#include <set>
#include <armadillo>



// Forward declaration of InteractorStyle
class InteractorStyle;

// Forward declaration of Mainwindow
class Mainwindow;


/**
Enum defining the state "transform direction" drop down list
*/
enum class TransformDirection {RADIAL, NORMAL_POINT,NORMAL_AVERAGED};


/**
Enum defining the state "interpolation type" drop down list
*/
enum class InterpolationType {UNIFORM, LINEAR, PARABOLIC};


/**
Enum defining the state "transform selection" drop down list
*/
enum class TransformSelection {SELECTED, NCLOSEST};




/**
Declaration of the ModifyAreaWidget class. ModifyAreaWidget refers to the
widget displayed on the screen when the user selects at least one vertex of the
displayed shape model by means of the rectangular box selector. The widget that is then
displayed lists the IDs of the selected vertices, as well as a choice of possible geometric
transforms to be applied to them. A horizontal slider_magnitude enalbes the user to choose
the magnitude of the transform
*/
class ModifyAreaWidget : public QDialog {
	Q_OBJECT

public:
	/**
	Constructor.
	@param parent Pointer to parent widget (here, pointer to instance of Mainwindow )
	@param interactor_style Pointer to the interactor accessing the shape model currently displayed
	*/

	ModifyAreaWidget(Mainwindow * parent,
	InteractorStyle * interactor_style);

	/**
	Data Setter The pointer passed as arguments allow the widget to have access to the
	point properties
	@param interactor_style Pointer to the interactor accessing the shape model currently displayed
	*/
	void set_data(InteractorStyle * interactor_style);

	/**
	Highlights the cells corresponding to the selected points
	*/
	void highlight_selected_cells();

	/**
	Computes the normals of the selected cells
	*/
	void compute_selected_cells_normals();

	/**
	Finds the position of the average of the selected vertices and 
	the coordinates of the vertex that is closest to this average
	*/
	void find_blob_center();


	void close();

	QTableWidget * table;

	QHBoxLayout * slider_magnitude_layout;
	QHBoxLayout * slider_neighbors_layout;

	QVBoxLayout * main_layout;

	QSlider * slider_magnitude;
	QSlider * slider_neighbors;

	QLineEdit * slider_magnitude_value;
	QLineEdit * slider_neighbors_value;

	QDialogButtonBox * button_box;

	QWidget * slider_magnitude_holder_widget;
	QWidget * slider_neighbors_holder_widget;

	QLabel * transform_direction_title;
	QLabel * interpolation_type_title;
	QLabel * transform_selection_title;

	QLabel * slider_magnitude_title;
	QLabel * slider_neighbors_title;

	QComboBox * transform_direction_list;
	QComboBox * interpolation_type_list;
	QComboBox * transform_selection_list;

	// Slots
private slots:
	/**
	Displays the current slider_magnitude position in the corresponding QEditText widget
	@param pos Slider position
	*/
	void show_new_slider_magnitude_pos(int pos);

	/**
	Displays the current slider_neighbors position in the corresponding QEditText widget
	@param pos Slider position
	*/
	void show_new_slider_neighbors_pos(int pos);


	/**
	Redraw the shape model based on the slider_magnitude position
	@param pos Slider position
	*/
	void update_view(int pos);

	/**
	Conveniency slot called when the view must be updated 
	without changing the slider_magnitude position
	*/
	void update_view_unchanged_slider_magnitude();

	/**
	Finds the indices of the N closest neighboring vertices to the selected blob center
	and uses this piece of information to update the view accordingly (using the 
	N-closest neighbor selection)
	@param N Number of neighboring vertices to identify
	*/
	void find_N_neighbors_indices_and_update_view(const int N);

	/**
	Sets the position slider_magnitude value via the QLineEdit field
	*/
	void set_new_slider_magnitude_pos();

	/**
	Sets the position slider_neighbors value via the QLineEdit field
	*/
	void set_new_slider_neighbors_pos();

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

	/**
	Finds the indice of the N closest neighbors to the blob center
	@param N Number of neighbors to identify
	*/
	void find_N_neighbors_indices(const int N);

	QStringList labels;

	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	vtkSmartPointer<vtkPolyData> active_selected_points_polydata;

	vtkSmartPointer<vtkPolyData> all_points_polydata;
	vtkSmartPointer<vtkPolyData> selected_cells_polydata;
	vtkSmartPointer<vtkPolyData> unselected_cells_polydata;

	vtkSmartPointer<vtkPoints> selected_points;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	
	vtkSmartPointer<vtkIdTypeArray> selected_polys_ids;
	vtkSmartPointer<vtkIdTypeArray> unselected_polys_ids;
	vtkSmartPointer<vtkPoints> new_selected_points_coordinates;

	vtkSmartPointer<vtkIdTypeArray> polys_ids;
	vtkSmartPointer<vtkIdTypeArray> visible_points_global_ids_from_local_index;

	vtkSmartPointer<vtkIdList> cell_ids;
	vtkDataArray * selected_cells_normals;
	vtkSmartPointer<vtkDoubleArray > averaged_normal_array;

	TransformDirection transform_direction;
	InterpolationType interpolation_type;
	TransformSelection transform_selection;

	arma::vec blob_center_position;
	arma::vec blob_average_position;
	int blob_center_id;

	vtkSmartPointer<vtkIdList> N_closest_vertices_indices;
	vtkSmartPointer<vtkKdTreePointLocator> point_locator;

	Mainwindow * parent;
};




#endif

