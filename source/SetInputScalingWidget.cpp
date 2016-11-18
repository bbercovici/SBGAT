#include "SetInputScalingWidget.hpp"

SetInputScalingWidget::SetInputScalingWidget(Mainwindow * parent) {
	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;

	this -> main_layout = new QVBoxLayout();
	this -> wrapping_layout = new QGridLayout();

	this -> scaling_factor_qlineedit = new QLineEdit(this);
	this -> input_unit_title = new QLabel("Input unit: ", this);
	this -> scaling_factor_title = new QLabel("Scaling factor: ", this);
	this -> input_unit_combobox = new QComboBox(this);
	this -> wrapping_widget = new QGroupBox("Scaling factor", this);
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, this);

	this -> input_unit_combobox -> addItem("meters (m)");
	this -> input_unit_combobox -> addItem("kilometers (km)");

	this -> wrapping_widget -> setLayout(this -> wrapping_layout);

	this -> wrapping_layout -> addWidget(input_unit_title, 0, 0, 1, 1);
	this -> wrapping_layout -> addWidget(input_unit_combobox, 0, 1, 1, 1);

	this -> wrapping_layout -> addWidget(scaling_factor_title, 1, 0, 1, 1);
	this -> wrapping_layout -> addWidget(scaling_factor_qlineedit, 1, 1, 1, 1);

	this -> setLayout(this -> main_layout);
	this -> main_layout -> addWidget(this -> wrapping_widget);
	this -> main_layout -> addWidget(this -> button_box);

	this -> main_layout -> addStretch(1);
	this -> setLayout(this -> main_layout);

	scaling_factor_qlineedit -> setText(QString::number(1));

	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(close()));
	connect(this -> input_unit_combobox, SIGNAL(currentIndexChanged(int)), this,
	        SLOT(update_scaling_factor()));

}


void SetInputScalingWidget::close() {
	this -> parent -> set_scaling_factor(scaling_factor_qlineedit -> text().toDouble());

	this -> parent -> lateral_dockwidget -> hide();

	QDialog::close();
}

void SetInputScalingWidget::update_scaling_factor() {
	switch (this -> input_unit_combobox -> currentIndex() ) {
	case 0:
		scaling_factor_qlineedit -> setText(QString::number(1));
		break;
	case 1:
		scaling_factor_qlineedit -> setText(QString::number(1000));
		break;

	}

}