#include <QApplication>
#include "Mainwindow.hpp"
#include <vtkVersion.h>

#include <ShapeModelImporter.hpp>


int main( int argc, char** argv ) {
	QSurfaceFormat surfaceFormat = QVTKOpenGLWidget::defaultFormat();

	surfaceFormat.setSamples(8);

	surfaceFormat.setStencilBufferSize(8);
	QSurfaceFormat::setDefaultFormat(surfaceFormat);
	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	Mainwindow window;

	return app.exec();


}
