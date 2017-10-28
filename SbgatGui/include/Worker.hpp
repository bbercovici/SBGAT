/*!
@file Worker.hpp
\author Benjamin Bercovici
\date July 22, 2017
\brief Stores definition of the Worker class.
*/

#ifndef HEADER_WORKER
#define HEADER_WORKER

#include <QObject>
#include <memory>
#include <DynamicAnalyses.hpp>
#include <chrono>
#include "ModelDataWrapper.hpp"

namespace SBGAT_GUI {



/*!
@class Worker
\author Benjamin Bercovici
\date July 22, 2017
\brief Worker class. Provides interface between SbgatCore, SbgatGUI and QThreads

\details This class is used an interface to a QThread, enabling heavy data processing outside of the
main thread within which the GUI lives. This allows the GUI to remain responsive. Its implementation
was based on the very informative blog post by Maya Posch at
https://mayaposch.wordpress.com/2011/11/01/how-to-really-truly-use-qthreads-the-full-explanation/
*/
class Worker: public QObject {
	Q_OBJECT

public:
	/**
	@param dyn_analyses pointer to DynamicAnalyses object instantiated from the shape model of interest
	@param mu standard gravitational parameter (m^3/s^-2)
	@param model_data pointer to the ModelDataWrapper representative of the shape being operated on
	@param name name of the shape
	*/
	Worker(std::shared_ptr<SBGAT_CORE::DynamicAnalyses> dyn_analyses,
	       double mu,
	       std::shared_ptr<ModelDataWrapper> model_data,
	       std::string name) ;

	virtual ~Worker() ;

public slots:
	/**
	Computes the global PGM acceleration over the provided shape model
	*/
	void process_pgm_acc();

	/**
	Computes the global PGM potentials over the provided shape model
	*/
	void process_pgm_pot() ;


signals:
	void progress(int);
	void finished();
	void error(QString err);
	/**
	Passes on a QString message to the console on the mainwindow
	@param log_out Log out messa
	*/
	void logging_out(QString log_out);

	/**
	Enables/disables log console
	@param console_status log console status
	*/
	void free_shape_table(bool console_status);

	/**
	Enables/disables menu bar
	@param menu_status menu bar status
	*/
	void free_menu_bar(bool menu_status);



private:

	std::shared_ptr<SBGAT_CORE::DynamicAnalyses> dyn_analyses;
	double mu;
	std::string name;
	std::shared_ptr<ModelDataWrapper> model_data;


};

}

#endif