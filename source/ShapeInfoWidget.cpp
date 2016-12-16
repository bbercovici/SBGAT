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
	
	unsigned int facets = this -> parent -> get_asteroid() -> get_polydata() -> GetNumberOfPolys();
	unsigned int vertices = this -> parent -> get_asteroid() -> get_polydata() -> GetNumberOfPoints();
	unsigned int edges = (unsigned int)(1.5 * (double)(facets));
	double length = this -> parent -> get_asteroid() -> get_polydata() -> GetLength();

	vtkSmartPointer<vtkMassProperties> mass_properties = vtkMassProperties::New();
	mass_properties -> SetInputData(this -> parent -> get_asteroid() -> get_polydata());
	mass_properties -> Update();

	std::string message =  "Vertices: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(vertices) ));

	message = "Facets: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(facets) ));

	message =  "Edges: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message + std::to_string(edges) ));

	this -> text_area -> appendPlainText(" ");

	message = "Characteristic length: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message) + QString::number(length,'e',5) + " m");

	message = "Area: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message)
	                                     + QString::number(mass_properties -> GetSurfaceArea ()/1e6,'e',5) + " km^2");

	message = "Volume: ";
	this -> text_area -> appendPlainText(QString::fromStdString(message)
	                                     + QString::number(mass_properties -> GetVolume ()/1e9,'e',5) + " km^3");



}

void ShapeInfoWidget::close() {


	this -> parent -> lateral_dockwidget -> hide();
	this -> parent -> set_action_status(true, this -> parent -> open_ComputePGMWidget_action);


	QDialog::close();
}