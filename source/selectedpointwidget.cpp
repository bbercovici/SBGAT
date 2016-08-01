#include "selectedpointwidget.h"

SelectedPointWidget::SelectedPointWidget(vtkSmartPointer<vtkPolyData> points_polydata,
        vtkSmartPointer<vtkPolyData> selected_points_polydata) {

	// The different GUI elements are created
	table = new QTableWidget();
	table -> setShowGrid(true);
	layout = new QHBoxLayout();
	list_holder_layout = new QVBoxLayout(this);
	list_holder_widget = new QWidget();
	transform_direction_list = new QComboBox();
	interpolation_type_list = new QComboBox();
	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Vertical);
	transform_direction_title = new QLabel("Transform direction:");
	interpolation_type_title = new QLabel("Interpolation type:");

	// The drop-down lists are added to the corresponding layout
	list_holder_widget -> setLayout(list_holder_layout);
	list_holder_layout -> setSpacing(0);
	list_holder_layout -> setMargin(0);

	list_holder_layout -> addWidget(transform_direction_title, Qt::AlignTop);
	list_holder_layout -> addWidget(transform_direction_list, Qt::AlignTop);
	list_holder_layout -> addWidget(interpolation_type_title, Qt::AlignTop);
	list_holder_layout -> addWidget(interpolation_type_list, Qt::AlignTop);





	// The two drop-down lists are filled
	transform_direction_list -> insertItem(0, "Radial");
	transform_direction_list -> insertItem(1, "Along blob normal");

	interpolation_type_list -> insertItem(0, "Uniform (0th order)");
	interpolation_type_list -> insertItem(1, "Linear (1st order)");
	interpolation_type_list -> insertItem(2, "Parabolic (2nd order)");

	// The selected points and the full point cloud are made accessible to the widget
	this -> selected_points_polydata = selected_points_polydata;
	this -> points_polydata = points_polydata;

	//  The dimensions of the table are set
	table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	table -> setColumnCount(4);

	labels << "ID" << "x" << "y" << "z";
	table -> setHorizontalHeaderLabels(labels);
	layout -> addWidget(table);
	layout -> addWidget(list_holder_widget);
	layout -> addWidget(button_box);

	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));

	this -> setLayout(layout);
	this -> populate();

}

void SelectedPointWidget::populate() {

	// Ids of selected points
	vtkSmartPointer<vtkDataArray> ids = this -> selected_points_polydata -> GetPointData() -> GetArray("ids");
	for (int row = 0; row < this -> selected_points_polydata -> GetNumberOfPoints(); ++row) {

		// The ids of each selected point is added in the first column
		QTableWidgetItem * id_item = new QTableWidgetItem(tr("%1").arg(* (ids -> GetTuple (row))));

		table -> setItem(row, 0, id_item);

		// The coordinates of each selected point are added to the other columns
		double p[3];
		this -> selected_points_polydata -> GetPoint(row, p);

		QTableWidgetItem * x_item = new QTableWidgetItem(tr("%1").arg(p[0]));
		QTableWidgetItem * y_item = new QTableWidgetItem(tr("%1").arg(p[1]));
		QTableWidgetItem * z_item = new QTableWidgetItem(tr("%1").arg(p[2]));

		table -> setItem(row, 1, x_item);
		table -> setItem(row, 2, y_item);
		table -> setItem(row, 3, z_item);

	}


}
