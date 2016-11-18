#include <QApplication>
#include "Mainwindow.hpp"

#include "Asteroid.hpp"
#include "Vect.hpp"
#include <fstream>

int main( int argc, char** argv ) {
	QApplication app( argc, argv );
	QLocale::setDefault(QLocale::c());


	std::ifstream file ("../resources/KW4alpha.grav");
	Asteroid asteroid = Asteroid(file);
	
	Vect Xc(3);
	Xc[0] = 1000;
	Xc[1] = 100;
	Xc[2] = 100;

	Vect grav = asteroid.PolyGrav(Xc, false);
	std::cout << grav[0] << " " << grav[1] << " " << grav[2] << std::endl;

	Mainwindow window;

	return app.exec();
}
