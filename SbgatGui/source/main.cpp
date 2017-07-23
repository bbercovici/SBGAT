#include <QApplication>
#include "Mainwindow.hpp"
#include <vtkVersion.h>

int main( int argc, char** argv ) {


	QSurfaceFormat surfaceFormat = QVTKOpenGLWidget::defaultFormat();

	QSurfaceFormat::setDefaultFormat(surfaceFormat);

	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	SBGAT_GUI::Mainwindow window;

	return app.exec();


}
