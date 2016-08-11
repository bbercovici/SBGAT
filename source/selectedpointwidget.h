
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
#include "vtkObjectFactory.h"


// Forward declaration of InteractorStyle
class InteractorStyle;

// Forward declaration of MainWindow
class MainWindow;



class vtkPolyData_tracked : public vtkPolyData{
public:
    static vtkPolyData_tracked * New();
    vtkTypeMacro(vtkPolyData_tracked, vtkPolyData);
    vtkPolyData_tracked();
    virtual ~vtkPolyData_tracked();
};


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
	Reset

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
	// vtkSmartPointer<vtkPolyData> selected_cells_polydata;
	// vtkSmartPointer<vtkPolyData> unselected_cells_polydata;


	// vtkSmartPointer<vtkPolyData_tracked> selected_points_polydata;
	// vtkSmartPointer<vtkPolyData_tracked> points_polydata;
	vtkSmartPointer<vtkPolyData_tracked> selected_cells_polydata;
	vtkSmartPointer<vtkPolyData_tracked> unselected_cells_polydata;


	vtkSmartPointer<vtkPoints> selected_points;
	std::vector<vtkSmartPointer<vtkActor> > actor_vector;
	vtkSmartPointer<vtkCellArray> selected_polys_cell_array;
	vtkSmartPointer<vtkCellArray> unselected_polys_cell_array;
	vtkSmartPointer<vtkIdTypeArray> selected_polys_ids;
	vtkSmartPointer<vtkIdTypeArray> unselected_polys_ids;
	vtkSmartPointer<vtkPoints> new_selected_points_coordinates;

	vtkSmartPointer<vtkIdTypeArray> polys_ids;
	vtkSmartPointer<vtkIdTypeArray> visible_v_ids;
	vtkSmartPointer<vtkIdList> cell_ids;


	MainWindow * mainwindow;
};




#endif

