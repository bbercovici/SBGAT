
#ifndef HEADER_COMPUTEPGMWIDGET
#define HEADER_COMPUTEPGMWIDGET

#include "Mainwindow.hpp"
#include "Asteroid.hpp"
#include "InteractorStyle.hpp"

#include <QVBoxLayout>
#include <QPushButton>
#include <QGridLayout>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QLabel>
#include <QLineEdit>
#include <QDoubleValidator>
#include <QPlainTextEdit>
#include <QErrorMessage>
#include <QComboBox>


#include <vtkModifiedBSPTree.h>
#include <vtkSphereSource.h>
#include <vtkActorCollection.h>
#include <vtkActor2DCollection.h>
#include <vtkLookupTable.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>


#include <armadillo>

// Forward declaration of Mainwindow
class Mainwindow;

class ComputePGMWidget : public QDialog {
	Q_OBJECT

public:

	/**
	Constructor.
	@param parent Pointer to parent widget (here, pointer to instance of Mainwindow )
	*/
	ComputePGMWidget(Mainwindow * parent);

	QVBoxLayout * main_layout;
	QGridLayout * physical_properties_layout;
	QGridLayout * spin_properties_box_layout;

	QGridLayout * compute_acceleration_layout;
	QGridLayout * compute_pgm_layout;
	QGridLayout * visualize_pgm_layout;

	QLabel * density_title_label;
	QLabel * density_unit_label;
	QLineEdit * density_qlineedit;


	QLabel * spin_rate_title_label;
	QLabel * spin_rate_unit_label;
	QLineEdit * spin_rate_qlineedit;

	QLabel * spin_axis_title_label;

	QLabel * spin_x_coordinate_label;
	QLabel * spin_y_coordinate_label;
	QLabel * spin_z_coordinate_label;

	QLineEdit * spin_x_coordinate_qlineedit;
	QLineEdit * spin_y_coordinate_qlineedit;
	QLineEdit * spin_z_coordinate_qlineedit;


	QLabel * x_coordinate_label;
	QLabel * y_coordinate_label;
	QLabel * z_coordinate_label;

	QLineEdit * x_coordinate_qlineedit;
	QLineEdit * y_coordinate_qlineedit;
	QLineEdit * z_coordinate_qlineedit;

	QPlainTextEdit * compute_acceleration_plainedit;

	QDialogButtonBox * button_box;
	QPushButton * compute_local_acceleration_button;
	QPushButton * compute_global_acceleration_button;
	QPushButton * load_global_acceleration_button;
	QPushButton * save_global_acceleration_button;

	QComboBox * visualization_combobox;

	QPushButton * show_visualization_button;

	QGroupBox * physical_properties_box;
	QGroupBox * compute_acceleration_box;
	QGroupBox * spin_properties_box;

	QGroupBox * compute_pgm_box;
	QGroupBox * visualize_pgm_box;

public slots:
	void close();

private slots:
	void compute_local_pgm();
	void compute_global_pgm();
	void load_global_pgm();
	void save_global_pgm();
	void move_local_acceleration_sphere();

	void clicked_visualization_button();

	void update_asteroid_state();


protected:
	Mainwindow * parent;
	vtkSmartPointer<vtkSphereSource> specific_point_sphere;
	vtkSmartPointer<vtkActor> specific_point_actor;
	vtkSmartPointer<vtkModifiedBSPTree> shape_mod_obb_tree;

	void cleanup();
	void show_acceleration_magnitude();
	void show_normal_acceleration_angle();
	void show_radial_acceleration_angle();
	void show_radial_acceleration_component();
	void show_normal_acceleration_component();
	void show_orthoradial_acceleration_magnitude();
	void show_orthonormal_acceleration_magnitude();
	void show_slopes();

};




#endif