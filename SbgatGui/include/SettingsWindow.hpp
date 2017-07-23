
#ifndef HEADER_SETTINGSWINDOW
#define HEADER_SETTINGSWINDOW

#include <QMainWindow>
#include <QVBoxLayout>
#include <QColorDialog>
#include <QColor>
#include <QGroupBox>
#include <QDialogButtonBox>
#include <QLabel>
#include <QComboBox>
#include <QCheckBox>

#include <QVTKOpenGLWidget.h>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {
	
class Mainwindow;

class SettingsWindow : public QDialog {
	Q_OBJECT

public:

	/**
	Constructor
	*/
	SettingsWindow(Mainwindow * parent) ;


private slots:

	void accept();
	void open_color_dialog();


protected:

	Mainwindow * parent;

	QComboBox * aa_frames_combo_box;
	QPushButton * open_color_dialog_button;
	QCheckBox * use_gradient_checkbox;


	double rgb_button[3];


};
}
#endif