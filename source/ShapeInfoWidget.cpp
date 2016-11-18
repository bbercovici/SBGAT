#include "ShapeInfoWidget.hpp"

ShapeInfoWidget::ShapeInfoWidget(Mainwindow * parent) {
	this -> setAttribute(Qt::WA_DeleteOnClose);
	this -> parent = parent;

	this -> main_layout = new QVBoxLayout();
	this -> text_area = new QPlainTextEdit(this);
	this -> text_area -> setReadOnly(true);
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok, Qt::Horizontal, this);

	this -> main_layout -> addWidget(this -> text_area);
	this -> main_layout -> addWidget(this -> button_box);

	this -> main_layout -> addStretch(1);
	this -> setLayout(this -> main_layout);

	this -> setupUI();

	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(close()));

}

void ShapeInfoWidget::setupUI() {
	InteractorStyle * mainwindow_interactor = static_cast< InteractorStyle * > (this -> parent -> get_render_window_interactor()
	        -> GetInteractorStyle());
	vtkPolyData * all_points_polydata = mainwindow_interactor -> get_all_points_polydata();
	unsigned int facets = all_points_polydata -> GetNumberOfPolys();
	unsigned int vertices = all_points_polydata -> GetNumberOfPoints();
	unsigned int edges = (unsigned int)(1.5 * (double)(facets));
	double length = all_points_polydata -> GetLength();
	std::string message =  "Vertices: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(vertices) ));

	message = "Facets: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(facets) ));

	message =  "Edges: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(edges) ));

	message = "Characteristic length: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(length) + " m"));

}

void ShapeInfoWidget::close() {


	this -> parent -> lateral_dockwidget -> hide();
	this -> parent -> set_action_status(true, this -> parent -> openComputePGMWidgetAct);


	QDialog::close();
}