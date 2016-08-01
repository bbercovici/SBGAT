#ifndef HEADER_NEWSHAPEMODELDIALOG
#define HEADER_NEWSHAPEMODELDIALOG

#include <QDialog>
#include <QHBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QIntValidator>
#include <QDialogButtonBox>
#include <iostream>

/**
Declaration fo the NewShapeModelDialog class. QDialog widget 
allowing the user to choose the number of normally distributed points in space 
used to create a convex shape model (for testing purposes).
*/


class NewShapeModelDialog : public QDialog {
	Q_OBJECT

public:
	NewShapeModelDialog();
	int result;
	QLineEdit * input_box;

	// Slots
private slots:


private:
	QHBoxLayout * layout;
	QLabel * label;
	QDialogButtonBox * button_box;



};


#endif
