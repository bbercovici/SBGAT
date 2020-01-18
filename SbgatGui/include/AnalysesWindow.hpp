/** MIT License

Copyright (c) 2018 Benjamin Bercovici and Jay McMahon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/


#ifndef HEADER_ANALYSESWINDOW
#define HEADER_ANALYSESWINDOW

#include <QMainWindow>
#include <QGroupBox>
#include <QComboBox>
#include <QCheckBox>
#include <QVBoxLayout>
#include <QPushButton>
#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>

#include "Mainwindow.hpp"

namespace SBGAT_GUI {

	class Mainwindow;

/*!
@class AnalysesWindow
\author Benjamin Bercovici
\date July 22, 2017
\brief AnalysesWindow a window superclass

*/

	class AnalysesWindow : public QDialog {
		Q_OBJECT

	public:

	/**
	Creates the settings window
	@param parent pointer to parent window.
	*/
		AnalysesWindow(Mainwindow * parent);

		Mainwindow * parent;
		
		QComboBox * prop_combo_box;
		
	protected:

		QVBoxLayout * analyses_layout;



	};
}
#endif