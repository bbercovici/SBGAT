#include "Worker.hpp"

using namespace SBGAT_GUI;


Worker::Worker(std::shared_ptr<SBGAT_CORE::DynamicAnalyses> dyn_analyses,
               double density,
               std::shared_ptr<ModelDataWrapper> model_data,
               std::string name) {

	this -> dyn_analyses = dyn_analyses;
	this -> density = density;
	this -> name = name;
	this -> model_data = model_data;

}

Worker::~Worker() {
}

void Worker::process_pgm_acc() {

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	this -> dyn_analyses -> compute_pgm_accelerations(density);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	model_data -> set_global_pgm_acc(true);
	model_data -> set_grav_slopes(false);


	// Log out
	QString log_out = QString::fromStdString("\n- Done computing in ")
	                  + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s");

	emit logging_out(log_out);
	emit free_shape_table(true);
	emit free_menu_bar(true);


	emit finished();
}


void Worker::process_pgm_pot() {

	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	this -> dyn_analyses -> compute_pgm_potentials(density);
	end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end - start;

	model_data -> set_global_pgm_pot(true);

	// Log out
	QString log_out = QString::fromStdString("\n- Done computing in ")
	                  + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s");

	emit logging_out(log_out);
	emit free_shape_table(true);
	emit free_menu_bar(true);

	emit finished();
}




