#include "newShapeModelDialog.h"

NewShapeModelDialog::NewShapeModelDialog() {
	layout = new QHBoxLayout();
	label = new QLabel("Number of points in generating point cloud");

	input_box = new QLineEdit();
	input_box -> setValidator( new QIntValidator(4, int(1e5), this) );

	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);

	layout -> addWidget(label);
	layout -> addWidget(input_box);
	layout -> addWidget(button_box);

	connect(button_box,SIGNAL(accepted()),this,SLOT(accept()));
	connect(button_box,SIGNAL(rejected()),this,SLOT(reject()));

	this -> setLayout(layout);
}


   




