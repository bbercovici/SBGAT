#include "selectedpointwidget.h"

SelectedPointWidget::SelectedPointWidget(vtkSmartPointer<vtkPolyData> points_polydata,
        vtkSmartPointer<vtkPolyData> selected_points_polydata,
        bool * widget_is_open, vtkSmartPointer<InteractorStyle> interactor_style) {

	// The different GUI elements are created

	layout = new QHBoxLayout();
	list_holder_layout = new QVBoxLayout();
	list_holder_widget = new QWidget();
	transform_direction_list = new QComboBox();
	interpolation_type_list = new QComboBox();
	button_box = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, Qt::Vertical);
	transform_direction_title = new QLabel("Transform direction:");
	interpolation_type_title = new QLabel("Interpolation type:");
	button_show_vertex_table = new QPushButton("Show vertex table");

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

	// The different widgets are added to the outer layout
	layout -> addWidget(list_holder_widget);
	layout -> addWidget(button_box, Qt::AlignCenter);
	layout -> addWidget(button_show_vertex_table, Qt::AlignCenter);
	this -> setLayout(layout);

	// The different buttons are connected to the corresponding slots
	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(reject()));
	connect(button_show_vertex_table, SIGNAL(clicked()), this, SLOT(show_vertex_table()));


	// This forces the widget to stay on top of all others
	Qt::WindowFlags flags = windowFlags();
	this -> setWindowFlags(flags | Qt::WindowStaysOnTopHint);

	// This prevents another instance of the widget to be opened
	this -> widget_is_open = widget_is_open;
	*this ->widget_is_open = true;

	// This links the selection widget to the interactor
	this -> interactor_style = interactor_style;

}

void SelectedPointWidget::show_vertex_table() {
	this -> table = new QTableWidget();
	this -> table -> setShowGrid(true);
	this -> table -> setEditTriggers(QAbstractItemView::NoEditTriggers); //  prevents the user from editing the items

	//  The dimensions of the table are set
	this -> table -> setRowCount(this -> selected_points_polydata -> GetNumberOfPoints());
	this -> table -> setColumnCount(4);

	this -> labels << "ID" << "x" << "y" << "z";
	this -> table -> setHorizontalHeaderLabels(labels);
	this -> populate_vertex_table();
	this -> table -> show();

}

void SelectedPointWidget::populate_vertex_table() {

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

void SelectedPointWidget::accept() {
	(*this -> widget_is_open) = false;
	this -> interactor_style -> transform_points();
	QDialog::accept();


}

void SelectedPointWidget::reject() {

	(*this -> widget_is_open) = false;
	QDialog::reject();


}
