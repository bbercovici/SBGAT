
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

#include <vtkModifiedBSPTree.h>
#include <vtkSphereSource.h>


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
	QGridLayout * compute_acceleration_layout;

	QLabel * density_title_label;
	QLabel * density_unit_label;
	QLineEdit * density_qlineedit;

	QLabel * gravitational_constant_title_label;
	QLabel * gravitational_constant_unit_label;
	QLineEdit * gravitational_constant_qlineedit;

	QLabel * scaling_factor_title_label;
	QLineEdit * scaling_factor_qlineedit;

	QLabel * x_coordinate_label;
	QLabel * y_coordinate_label;
	QLabel * z_coordinate_label;

	QLineEdit * x_coordinate_qlineedit;
	QLineEdit * y_coordinate_qlineedit;
	QLineEdit * z_coordinate_qlineedit;

	QPlainTextEdit * compute_acceleration_plainedit;

	QDialogButtonBox * button_box;
	QPushButton * compute_PGM_button;
	QGroupBox * physical_properties_box;
	QGroupBox * compute_acceleration_box;

public slots:
	void close();

private slots:
	void compute_pgm();
	void move_specific_point_sphere();


protected:
	Mainwindow * parent;
	vtkSmartPointer<vtkSphereSource> specific_point_sphere;
	vtkSmartPointer<vtkActor> specific_point_actor;

};




#endif