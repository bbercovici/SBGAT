// #include <QApplication>
// #include "Mainwindow.hpp"
// #include <vtkVersion.h>

#include "ShapeModel.hpp"
#include "DynamicAnalyses.hpp"

int main( int argc, char** argv ) {

	// QApplication app( argc, argv );
	// QLocale::setDefault(QLocale::c());

	// Mainwindow window;

	// return app.exec();

	ShapeModel shape_model;
	shape_model.load("../resources/KW4alpha.obj");
	shape_model.check_normals_consistency();


	DynamicAnalyses dynamic_analyses(&shape_model);
	dynamic_analyses.compute_pgm();


	return 0;
}
