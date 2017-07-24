#ifndef HEADER_WORKER
#define HEADER_WORKER

#include <QObject>
#include <QPlainTextEdit>
#include <memory>
#include <DynamicAnalyses.hpp>
#include <chrono>
#include "ModelDataWrapper.hpp"

namespace SBGAT_GUI {



class Worker: public QObject {
	Q_OBJECT

public:
	Worker(std::shared_ptr<SBGAT_CORE::DynamicAnalyses> dyn_analyses,
	       double density,
	       std::shared_ptr<ModelDataWrapper> model_data,
	       std::string name,
	       QPlainTextEdit * log_console) ;
	virtual ~Worker() ;

public slots:
	void process_pgm_acc();
	void process_pgm_pot() ;


signals:
	void progress(int);
	void finished();
	void error(QString err);
	void logging_out(QString log_out);
	void free_gui(bool gui_status);



private:

	std::shared_ptr<SBGAT_CORE::DynamicAnalyses> dyn_analyses;
	double density;
	std::string name;
	QPlainTextEdit * log_console;
	std::shared_ptr<ModelDataWrapper> model_data;


};

}

#endif