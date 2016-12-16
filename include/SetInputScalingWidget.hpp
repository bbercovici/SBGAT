#ifndef HEADER_SETINPUTSCALINGWIDGET
#define HEADER_SETINPUTSCALINGWIDGET

#include "Mainwindow.hpp"
#include "Asteroid.hpp"
#include "InteractorStyle.hpp"

#include <QVBoxLayout>
#include <QPushButton>
#include <QGridLayout>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QLabel>
#include <QComboBox>

#include <armadillo>

// Forward declaration of Mainwindow
class Mainwindow;

class SetInputScalingWidget : public QDialog {
	Q_OBJECT

public:

	/**
	Constructor.
	@param parent Pointer to parent widget (here, pointer to instance of Mainwindow )
	*/
	SetInputScalingWidget(Mainwindow * parent);

public slots:
	void close();

private slots:
	/**
	Updates the scaling factor accordingly
	when a new unit is chosen
	*/
	void update_scaling_factor();

protected:
	Mainwindow * parent;

	QVBoxLayout * main_layout;
	QGridLayout * wrapping_layout;

	QComboBox * input_unit_combobox;
	QLineEdit * scaling_factor_qlineedit;
	QLabel * input_unit_title;
	QLabel * scaling_factor_title;
	QGroupBox * wrapping_widget;

	QDialogButtonBox * button_box;



};


#endif