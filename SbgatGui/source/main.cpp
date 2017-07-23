#include <QApplication>
#include "Mainwindow.hpp"
#include <vtkVersion.h>

#include <ShapeModelImporter.hpp>


int main( int argc, char** argv ) {


	QSurfaceFormat surfaceFormat = QVTKOpenGLWidget::defaultFormat();

	QSurfaceFormat::setDefaultFormat(surfaceFormat);

	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	Mainwindow window;

	return app.exec();


}
