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

#include "LCVisualizer.hpp"
#include <SBGATObsLightcurve.hpp>

#include <vtkGenericOpenGLRenderWindow.h>
#include <QHBoxLayout>
#include <QVBoxLayout>

#include <vtkActor2DCollection.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>
#include <vtkChartHistogram2D.h>
#include <vtkContextScene.h>
#include <vtkStringArray.h>
#include <vtkPlot.h>
#include <vtkColorTransferFunction.h>
#include <vtkAxis.h>
#include <vtkFloatArray.h>
#include <vtkTable.h>
#include <QLabel>


using namespace SBGAT_GUI;

LCVisualizer::LCVisualizer(LCWindow * parent,const std::vector<std::array<double, 2> > & measurements) {

	this -> parent = parent;
	this -> setWindowTitle("Visualize Light Curve");

	this -> qvtkWidget = new QVTKOpenGLWidget(this);
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);
	QVBoxLayout * visualizer_window_layout = new QVBoxLayout(this);
	visualizer_window_layout -> addWidget(this -> qvtkWidget,1);
	visualizer_window_layout -> addWidget(this -> button_box,0);

	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));


	this -> init(measurements);

	
}


void LCVisualizer::init(const std::vector<std::array<double, 2> > & measurements){

	// Creating view, histogram, transfer function and rendering window.
	vtkSmartPointer<vtkGenericOpenGLRenderWindow> render_window = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
	
	// vtkSmartPointer<vtkChartHistogram2D> histo = vtkSmartPointer<vtkChartHistogram2D>::New();
	// vtkSmartPointer<vtkColorTransferFunction> fun = vtkSmartPointer<vtkColorTransferFunction>::New();

	this -> view = vtkSmartPointer<vtkContextView>::New ();
	this -> view -> SetRenderWindow(render_window);
	this -> qvtkWidget -> SetRenderWindow(render_window);

	// Creating the chart
	vtkSmartPointer<vtkChartXY> chart =
	vtkSmartPointer<vtkChartXY>::New();
	chart->SetShowLegend(true);

	vtkSmartPointer<vtkTable> table =
	vtkSmartPointer<vtkTable>::New();

	vtkSmartPointer<vtkFloatArray> time =
	vtkSmartPointer<vtkFloatArray>::New();
	time->SetName("Time since epoch (days)");
	table->AddColumn(time);

	vtkSmartPointer<vtkFloatArray> luminosity = vtkSmartPointer<vtkFloatArray>::New();
	luminosity -> SetName("Luminosity");
	
	table -> AddColumn(luminosity);
	table -> SetNumberOfRows(measurements.size());

	// Extracting the max luminosity for normalization purposes
	double max_luminosity = measurements[0][1];

	for (int i = 0; i < measurements.size(); ++i){
		max_luminosity = std::max(max_luminosity,measurements[i][1]);
	}

	for (int i = 0; i < measurements.size(); ++i){
		table -> SetValue(i, 0, (measurements[i][0] - measurements[0][0]) / 86400);
		table -> SetValue(i, 1, measurements[i][1] / max_luminosity);

	}


	vtkPlot * points = chart -> AddPlot(vtkChart::POINTS);

	points -> SetInputData(table, 0, 1);

	chart -> GetAxis( vtkAxis::LEFT ) -> SetTitle("Luminosity");
	chart -> GetAxis( vtkAxis::BOTTOM ) -> SetTitle("Time since epoch (days)");
	chart -> GetAxis( vtkAxis::LEFT ) -> SetVisible(1);
	chart -> GetAxis( vtkAxis::BOTTOM ) -> SetVisible(1);

	this -> view -> GetScene() -> AddItem(chart);
		
	this -> resize(500, 500);

	this -> qvtkWidget -> update();
	this -> view -> GetRenderWindow() -> Render();
	this -> qvtkWidget -> GetRenderWindow() -> Render();
	this -> qvtkWidget -> repaint();
	this -> view -> GetInteractor()->Start();


}

