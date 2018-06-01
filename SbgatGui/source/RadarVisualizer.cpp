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

#include "RadarVisualizer.hpp"
#include <SBGATObsRadar.hpp>

#include <vtkGenericOpenGLRenderWindow.h>
#include <QHBoxLayout>
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

#include <QLabel>


using namespace SBGAT_GUI;

RadarVisualizer::RadarVisualizer(RadarWindow * parent,
	const std::vector<vtkSmartPointer<vtkImageData>> & images) {

	this -> parent = parent;
	this -> setWindowTitle("Visualize Radar Images");

	this -> qvtkWidget = new QVTKOpenGLWidget(this);
	// this -> qvtkWidget = new QVTKWidget(this);

	this -> images = images;
	this -> previous_image_button = new QPushButton("<", this);
	this -> next_image_button = new QPushButton(">", this);
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	QVBoxLayout * visualizer_window_layout = new QVBoxLayout(this);

	QWidget * navigation_widget = new QWidget(this);
	QHBoxLayout * navigation_widget_layout = new QHBoxLayout(navigation_widget);

	navigation_widget_layout -> addWidget(this -> previous_image_button);
	navigation_widget_layout -> addWidget(this -> next_image_button);
	navigation_widget_layout -> addWidget(this -> button_box);

	visualizer_window_layout -> addWidget(this -> qvtkWidget,1);
	visualizer_window_layout -> addWidget(navigation_widget,0);


	connect(this -> previous_image_button,SIGNAL(clicked()), this, SLOT(previous_image()));
	connect(this -> next_image_button,SIGNAL(clicked()), this, SLOT(next_image()));
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));

	this -> init();

	
}


void RadarVisualizer::init(){

	// Creating view, histogram, transfer function and rendering window.
	vtkSmartPointer<vtkGenericOpenGLRenderWindow> render_window = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
	
	vtkSmartPointer<vtkChartHistogram2D> histo = vtkSmartPointer<vtkChartHistogram2D>::New();
	vtkSmartPointer<vtkColorTransferFunction> fun = vtkSmartPointer<vtkColorTransferFunction>::New();

	this -> view = vtkSmartPointer<vtkContextView>::New ();
	this -> view -> SetRenderWindow(render_window);
	this -> qvtkWidget -> SetRenderWindow(render_window);

	// Creating the image
	vtkSmartPointer<vtkImageData> image = this -> images[this -> current_image_index];

	// Normalizing the image
	vtkDataArray * scalars = image -> GetPointData() -> GetScalars();
	double max_val = -1;

	for (vtkIdType tupleIdx = 0; tupleIdx < scalars -> GetNumberOfTuples(); ++tupleIdx){
		max_val = std::max(scalars -> GetTuple1(tupleIdx),max_val);
	}

	this -> resize(500, 500);

	// Setting the transfer function
	fun -> AddRGBPoint(0,0.0, 0.0, 0.0);
	fun -> AddRGBPoint(max_val,  1.0, 1.0, 1.0);
	fun -> Build();

	// Computing the histogram and setting axes titles
	histo -> SetInputData(image);
	histo -> SetTransferFunction(fun);


	histo -> GetAxis( vtkAxis::LEFT) -> SetTitle("Range bin");
	histo -> GetAxis( vtkAxis::BOTTOM) -> SetTitle("Range-rate bin");
	histo -> GetAxis( vtkAxis::LEFT) -> SetVisible(1);
	histo -> GetAxis( vtkAxis::BOTTOM) -> SetVisible(1);

	// std::cout << histo -> GetAxis( vtkAxis::LEFT) -> GetTitle() << std::endl;

    // this -> view -> GetScene() -> ClearItems();
	this -> view -> GetScene() -> AddItem(histo);
	
	// this -> qvtkWidget -> setMinimumSize(int((double)(max_x)/2),int((double)(max_y)/2));
	this -> qvtkWidget -> update();
	this -> view -> GetRenderWindow() -> Render();
	this -> qvtkWidget -> GetRenderWindow() -> Render();
	this -> qvtkWidget -> repaint();
	this -> view->GetInteractor()->Start();


}

void RadarVisualizer::next_image(){

	if (this -> current_image_index + 1 == this -> images.size()){
		this -> current_image_index = 0;
	}
	else{
		++this -> current_image_index;
	}

	auto image = this -> images[this -> current_image_index];
	auto histo = vtkChartHistogram2D::SafeDownCast(this -> view -> GetScene() -> GetItem(0));


	vtkDataArray * scalars = image -> GetPointData() -> GetScalars();
	double max_val = -1;

	for (vtkIdType tupleIdx = 0; tupleIdx < scalars -> GetNumberOfTuples(); ++tupleIdx){
		max_val = std::max(scalars -> GetTuple1(tupleIdx),max_val);
	}

	// Setting the transfer function
	vtkSmartPointer<vtkColorTransferFunction> fun = vtkSmartPointer<vtkColorTransferFunction>::New();

	fun -> AddRGBPoint(0,0.0, 0.0, 0.0);
	fun -> AddRGBPoint(max_val,  1.0, 1.0, 1.0);
	fun -> Build();

	// Computing the histogram and setting axes titles
	histo -> SetInputData(image);
	histo -> SetTransferFunction(fun);

	this -> qvtkWidget -> update();
	this -> qvtkWidget -> GetRenderWindow() -> Render();


}

void RadarVisualizer::previous_image(){


	if (this -> current_image_index == 0){
		this -> current_image_index = (int)(this -> images.size()) - 1;
	}
	else{
		--this -> current_image_index;
	}

	auto image = this -> images[this -> current_image_index];
	auto histo = vtkChartHistogram2D::SafeDownCast(this -> view -> GetScene() -> GetItem(0));

	vtkDataArray * scalars = image -> GetPointData() -> GetScalars();
	double max_val = -1;

	for (vtkIdType tupleIdx = 0; tupleIdx < scalars -> GetNumberOfTuples(); ++tupleIdx){
		max_val = std::max(scalars -> GetTuple1(tupleIdx),max_val);
	}

	// Setting the transfer function
	vtkSmartPointer<vtkColorTransferFunction> fun = vtkSmartPointer<vtkColorTransferFunction>::New();

	fun -> AddRGBPoint(0,0.0, 0.0, 0.0);
	fun -> AddRGBPoint(max_val,  1.0, 1.0, 1.0);
	fun -> Build();

	// Computing the histogram and setting axes titles
	histo -> SetInputData(image);
	histo -> SetTransferFunction(fun);

	this -> qvtkWidget -> update();
	this -> qvtkWidget -> GetRenderWindow() -> Render();




	
}

