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

#include <vtkDataArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>

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
	Constructor. The pointer passed as arguments allow the widget to have access to the
	point properties
	@param points_polydata Pointer to the vtkPolyData storing the vertices of the shape model
	@param selected_points_polydata Pointer to the vtkPolyData storing the selected vertices
	of the displayed shape model
	*/
	SelectedPointWidget(vtkSmartPointer<vtkPolyData> points_polydata,
	                    vtkSmartPointer<vtkPolyData> selected_points_polydata);


	QTableWidget * table;
	QHBoxLayout * layout;
	QVBoxLayout * list_holder_layout;

	QDialogButtonBox * button_box;
	QWidget * list_holder_widget;
	QLabel * transform_direction_title;
	QLabel * interpolation_type_title;

	QComboBox * transform_direction_list;
	QComboBox * interpolation_type_list;

	/**
	Populates the QTableWidget table with the relevant data
	*/
	void populate();

	// Slots
private slots:

private:
	void createActions();
	void createMenus();

	QStringList labels;
	vtkSmartPointer<vtkPolyData> selected_points_polydata;
	vtkSmartPointer<vtkPolyData> points_polydata;

};




#endif

