#include <QApplication>
#include "Mainwindow.hpp"
#include <vtkVersion.h>

#include <ShapeModelImporter.hpp>


int main( int argc, char** argv ) {

	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	Mainwindow window;

	return app.exec();


}
