#ifndef HEADER_SHAPEINFOWIDGET
#define HEADER_SHAPEINFOWIDGET

#include "Mainwindow.hpp"
#include "Asteroid.hpp"
#include "InteractorStyle.hpp"

#include <QVBoxLayout>
#include <QDialogButtonBox>
#include <QGroupBox>
#include <QPlainTextEdit>

#include <armadillo>

// Forward declaration of Mainwindow
class Mainwindow;

class ShapeInfoWidget : public QDialog {
	Q_OBJECT

public:

	/**
	Constructor.
	@param parent Pointer to parent widget (here, pointer to instance of Mainwindow )
	*/
	ShapeInfoWidget(Mainwindow * parent);

public slots:
	void close();

private slots:


protected:
	Mainwindow * parent;

	QVBoxLayout * main_layout;
	QPlainTextEdit * text_area;
	QDialogButtonBox * button_box;
	void setupUI();



};


#endif