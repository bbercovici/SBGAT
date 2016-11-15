#include <QApplication>
#include "Mainwindow.hpp"

int main( int argc, char** argv ) {
	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	Mainwindow window;

	return app.exec();
}
