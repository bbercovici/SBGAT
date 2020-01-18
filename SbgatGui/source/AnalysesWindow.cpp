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

#include "SHARMWindow.hpp"
#include <SBGATSphericalHarmo.hpp>

using namespace SBGAT_GUI;

AnalysesWindow::AnalysesWindow(Mainwindow * parent)  {

	this -> parent = parent;
	this -> analyses_layout = new QVBoxLayout(this);

	QGroupBox * shape_model_group = new QGroupBox(tr("Shape"));
	
	QGridLayout * shape_model_group_layout = new QGridLayout(shape_model_group);
	
	QLabel * shape_label = new QLabel("Shape model",this);

	this -> prop_combo_box = new QComboBox (this);

	shape_model_group_layout -> addWidget(shape_label,0,0,1,1);
	shape_model_group_layout -> addWidget(this -> prop_combo_box,0,1,1,1);

	this -> analyses_layout -> addWidget(shape_model_group);
	
}
