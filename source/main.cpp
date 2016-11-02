#include <QApplication>
#include "Mainwindow.hpp"

int main( int argc, char** argv ) {
	QApplication app( argc, argv );
	Mainwindow window;
	
	return app.exec();
}
