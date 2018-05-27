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
#include <vtkImageMapper.h>
#include <vtkActor2D.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <QHBoxLayout>
#include <vtkActor2DCollection.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkImageData.h>



#include <QLabel>



using namespace SBGAT_GUI;

RadarVisualizer::RadarVisualizer(RadarWindow * parent,
	const std::vector<vtkSmartPointer<vtkImageData>> & images) {

	this -> parent = parent;
	this -> setWindowTitle("Visualize Radar Images");

	this -> qvtkWidget = new QVTKOpenGLWidget(this);
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

	vtkSmartPointer<vtkImageMapper> imageMapper = vtkSmartPointer<vtkImageMapper>::New();
	imageMapper -> SetInputData(this -> images[0]);

	vtkSmartPointer<vtkActor2D> imageActor = vtkSmartPointer<vtkActor2D>::New();
	imageActor -> SetMapper(imageMapper);

   // Setup renderer
	this -> renderer = vtkSmartPointer<vtkRenderer>::New();

   // Setup render window
	vtkSmartPointer<vtkGenericOpenGLRenderWindow> render_window = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
	this -> qvtkWidget -> SetRenderWindow(render_window);
	this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(this -> renderer);

	render_window -> AddRenderer(renderer);
	this -> renderer -> AddActor2D(imageActor);

	int dimensions[3];
	int max_x = -1;
	int max_y = -1;

	for (int i = 0; i < this -> images.size(); ++i){
		this -> images[i] -> GetDimensions(dimensions);
		max_x = std::max(max_x,dimensions[0]);
		max_y = std::max(max_y,dimensions[1]);
	}

	this -> qvtkWidget -> setMinimumSize(int((double)(max_x)/2),int((double)(max_y)/2));
	this -> qvtkWidget -> update();
	this -> qvtkWidget -> repaint();
	this -> qvtkWidget -> GetRenderWindow() -> Render();

}

void RadarVisualizer::next_image(){

	vtkActor2DCollection * collection = this -> renderer -> GetActors2D();
	collection -> InitTraversal();

	vtkActor2D * actor =  collection -> GetNextActor2D();
	vtkImageMapper * mapper = vtkImageMapper::SafeDownCast(actor -> GetMapper());
	vtkSmartPointer<vtkImageData> current_image = mapper -> GetInput();

	auto current_image_iterator = std::find(this -> images.begin(), this -> images.end(), current_image);

	if (current_image_iterator != --this -> images.end()){
		mapper -> SetInputData(*(++current_image_iterator));	
	}
	else{
		mapper -> SetInputData(*this -> images.begin());	
	}


	mapper -> Update();
	
	this -> qvtkWidget -> update();
	this -> qvtkWidget -> GetRenderWindow() -> Render();

}

void RadarVisualizer::previous_image(){

	vtkActor2DCollection * collection = this -> renderer -> GetActors2D();
	collection -> InitTraversal();
	vtkActor2D * actor =  collection -> GetNextActor2D();

	vtkImageMapper * mapper = vtkImageMapper::SafeDownCast(actor -> GetMapper());
	vtkSmartPointer<vtkImageData> current_image = mapper -> GetInput();
	auto current_image_iterator = std::find(this -> images.begin(), this -> images.end(), current_image);

	if (current_image_iterator != this -> images.begin()){
		mapper -> SetInputData(*(--current_image_iterator));	
	}
	else{
		mapper -> SetInputData(*(--this -> images.end()));	
	}

	mapper -> Update();
	
	this -> qvtkWidget -> update();
	this -> qvtkWidget -> GetRenderWindow() -> Render();


	
}

