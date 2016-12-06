#include <QApplication>
#include "Mainwindow.hpp"

#include "Asteroid.hpp"
#include "Vect.hpp"
#include <fstream>

int main( int argc, char** argv ) {
	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());

	Mainwindow window;

	return app.exec();
}
