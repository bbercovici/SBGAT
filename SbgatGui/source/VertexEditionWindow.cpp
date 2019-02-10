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

#include "VertexEditionWindow.hpp"

#include <QVBoxLayout>
#include <QGroupBox>
#include <QPushButton>
#include <QCheckBox>
#include <QSpinBox>
#include <QDialogButtonBox>
#include <QTableWidget>
#include <QLabel>
#include <vtkCamera.h>


#include <vtkLightCollection.h>
#include <vtkLight.h>
#include <vtkShadowMapPass.h>

#include <vtkProperty.h>
#include <vtkLineSource.h>
#include <vtkSphereSource.h>

#include <vtkDoubleArray.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>


using namespace SBGAT_GUI;

VertexEditionWindow::VertexEditionWindow(Mainwindow * parent) : QDialog(parent,Qt::WindowStaysOnTopHint) {

	this -> parent = parent;
	this -> setWindowTitle("Vertex Edition");

	QVBoxLayout * main_layout = new QVBoxLayout(this);
	QGroupBox * slider_widget = new QGroupBox("Vertex motion",this);
	QGroupBox * neighborhood_widget = new QGroupBox("Neighborhood",this);

	QHBoxLayout * slider_layout = new QHBoxLayout(slider_widget);
	QHBoxLayout * neighborhood_layout = new QHBoxLayout(neighborhood_widget);
	this -> neighborhood_spin_box = new QSpinBox(this);

	this -> neighborhood_spin_box -> setMinimum(1);
	this -> neighborhood_spin_box -> setMaximum(999);

	this -> direction_combo_box = new QComboBox(this);
	this -> position_slider = new QSlider(Qt::Horizontal,this);

	slider_layout -> addWidget(this -> direction_combo_box);
	slider_layout -> addWidget(this -> position_slider);

	QLabel * neighborhood_size_label = new QLabel("Size",this);

	neighborhood_layout -> addWidget(neighborhood_size_label);
	neighborhood_layout -> addWidget(neighborhood_spin_box);


	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok| QDialogButtonBox::Cancel);
	main_layout -> addWidget(slider_widget);
	main_layout -> addWidget(neighborhood_widget);
	main_layout -> addWidget(button_box);
	this -> line_actor = vtkSmartPointer<vtkActor>::New();
	this -> line_mapper = vtkSmartPointer<vtkDataSetMapper>::New();

	this -> init();
	this -> update_direction(0);

	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this, SIGNAL(rejected()),this,SLOT(close()));
	connect(this -> position_slider, SIGNAL(sliderMoved(int)),this, SLOT(move_vertex()));
	connect(this -> direction_combo_box, SIGNAL(currentIndexChanged(int)),this, SLOT(update_direction(int)));
	connect(this -> neighborhood_spin_box,SIGNAL(valueChanged(int)),this,SLOT(update_neighborhood_size(int)));

}

void VertexEditionWindow::init(){

	this -> direction_combo_box -> insertItem(this -> direction_combo_box -> count(),QString::fromStdString("Normal"));
	this -> direction_combo_box -> insertItem(this -> direction_combo_box -> count(),QString::fromStdString("X"));
	this -> direction_combo_box -> insertItem(this -> direction_combo_box -> count(),QString::fromStdString("Y"));
	this -> direction_combo_box -> insertItem(this -> direction_combo_box -> count(),QString::fromStdString("Z"));

	PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetQueriedPoint(this -> queried_point);
	
	// this -> position_slider -> setMinimum(-10);
	// this -> position_slider -> setMaximum(-10);
	this -> position_slider -> setValue(50);

	this -> line_actor -> SetMapper(this -> line_mapper);
	this -> parent -> get_renderer() -> AddActor(this -> line_actor);



}

void VertexEditionWindow::update_direction(int new_index){








	switch(new_index){

		case 0:
		PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetNormalAtSelectedPoint(this -> dir);
		
		this -> line_actor -> GetProperty() -> SetColor(0.5,0.5,0.5);
		break;

		case 1:
		this -> dir[0] = 1;
		this -> dir[1] = 0;
		this -> dir[2] = 0;
		this -> line_actor -> GetProperty() -> SetColor(1,0,0);
		break;

		case 2:
		this -> dir[0] = 0;
		this -> dir[1] = 1;
		this -> dir[2] = 0;
		this -> line_actor -> GetProperty() -> SetColor(0,1,0);
		break;

		case 3:
		this -> dir[0] = 0;
		this -> dir[1] = 0;
		this -> dir[2] = 1;
		this -> line_actor -> GetProperty() -> SetColor(0,0,1);
		break;

	}

	double length = PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) ->  GetSelectedBodySize();

	double p0[3] = {
		this -> queried_point[0] + length * this -> dir[0], 
		this -> queried_point[1] + length * this -> dir[1], 
		this -> queried_point[2] + length * this -> dir[2]
	};
	double p1[3] = {
		this -> queried_point[0] - length * this -> dir[0], 
		this -> queried_point[1] - length * this -> dir[1], 
		this -> queried_point[2] - length * this -> dir[2]
	};

	vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
	lineSource -> SetPoint1(p0);
	lineSource -> SetPoint2(p1);
	lineSource -> Update();

	this -> line_mapper -> SetInputConnection(lineSource -> GetOutputPort());
	this -> line_actor -> GetProperty() -> SetLineWidth(4);
	this -> line_actor -> Modified();
	this -> move_vertex();
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

void VertexEditionWindow::move_vertex(){
	
	double length = PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) ->  GetSelectedBodySize();
	double val = 0.25 * length * (double(this -> position_slider -> value()) - 50)/50.;
	double translation_vector[3] = { 
		val * this -> dir[0], 
		val * this -> dir[1],  
		val * this -> dir[2]
	};

	double position[3] = {
		translation_vector[0],
		translation_vector[1],
		translation_vector[2]
	};
	PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetSelectedActor() -> SetPosition(position);
	PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetSelectedActor() -> Modified();
	this -> update_neighbors_positions();
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}

void VertexEditionWindow::close(){
	double zeros_pos[3] = {0,0,0};
	PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetSelectedActor() -> SetPosition(zeros_pos);
	PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetSelectedActor() -> Modified();
	
	for (int e = 0; e < this -> neighbors_actors_vector.size(); ++e){
		this -> parent -> get_renderer() -> RemoveActor(this -> neighbors_actors_vector[e]);
	}

	this -> neighbors_actors_vector.clear();
	this -> neighbors_indices_vector.clear();
	this -> neighbors_distances_vector.clear();

	this -> parent -> get_renderer() -> RemoveActor(this -> line_actor);
	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


	QDialog::close();
}

void VertexEditionWindow::accept(){

	int selected_row_index = this ->  parent -> prop_table -> selectionModel() -> currentIndex().row();
	std::string name = this -> parent ->  prop_table -> item(selected_row_index, 0) -> text() .toStdString();
	vtkSmartPointer<vtkPolyData> selected_shape = this -> parent -> get_wrapped_shape_data()[name]-> get_polydata();
	double length = PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) ->  GetSelectedBodySize();
	double val_0 = 0.25 * length * (double(this -> position_slider -> value()) - 50)/50.;
	
	double new_pos[3] = {
		this -> queried_point[0] + val_0 * this -> dir[0], 
		this -> queried_point[1] + val_0 * this -> dir[1], 
		this -> queried_point[2] + val_0 * this -> dir[2]
	};

	selected_shape -> GetPoints() -> SetPoint(PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) ->  GetQueriedPointId(),new_pos);

	for (int e = 0; e < this -> neighbors_distances_vector.size(); ++e ){
		double n[3];

		PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetNormalAtPoint(this -> neighbors_indices_vector[e],n);

		double val = std::exp(- std::pow(this -> neighbors_distances_vector[e]/(this -> neighbors_distances_vector.back()/3),2)) * val_0;
		double translation_vector[3] = { 
			val * n[0], 
			val * n[1],  
			val * n[2]
		};

		double old_pos[3];

		selected_shape -> GetPoints() -> GetPoint(this -> neighbors_indices_vector[e],old_pos);

		double new_pos[3] = {
			old_pos[0] + translation_vector[0],
			old_pos[1] + translation_vector[1],
			old_pos[2] + translation_vector[2]
		};

		selected_shape -> GetPoints() -> SetPoint(this -> neighbors_indices_vector[e],new_pos);

	}

	// Should update normals and kd tree here
	vtkSmartPointer<vtkPolyDataNormals> normals =
	vtkSmartPointer<vtkPolyDataNormals>::New();
	normals -> SetInputData(selected_shape);
	normals -> SplittingOff();
	normals -> ConsistencyOn();
	normals -> ComputeCellNormalsOff();
	normals -> ComputePointNormalsOn();
	normals -> Update();

	selected_shape -> GetPointData() -> GetArray("Normals") -> ShallowCopy(normals -> GetOutput() -> GetPointData() -> GetArray("Normals"));
	selected_shape -> GetPointData() -> GetArray("Normals") -> Modified();
	selected_shape -> GetPointData() -> Modified();
	selected_shape -> GetPoints() -> Modified();
	selected_shape -> Modified();


	this -> parent -> get_wrapped_shape_data()[name]-> get_tree() -> BuildLocator();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
    this -> parent -> prop_table -> setItem(selected_row_index, 1, new QTableWidgetItem("Modified"));

	this -> close();

}

void VertexEditionWindow::update_neighborhood_size(int neighborhood_size){

	vtkSmartPointer<vtkIdList> result = vtkSmartPointer<vtkIdList>::New();
	int selected_row_index = this ->  parent -> prop_table -> selectionModel() -> currentIndex().row();
	std::string name = this -> parent ->  prop_table -> item(selected_row_index, 0) -> text() .toStdString();
	this -> parent -> get_wrapped_shape_data()[name]-> get_tree() -> FindClosestNPoints(neighborhood_size, this -> queried_point, result);

	for (int e = 0; e < this -> neighbors_actors_vector.size(); ++e){
		this -> parent -> get_renderer() -> RemoveActor(this -> neighbors_actors_vector[e]);
	}

	this -> neighbors_actors_vector.clear();
	this -> neighbors_indices_vector.clear();
	this -> neighbors_distances_vector.clear();
	
	for (int e = 0; e < result -> GetNumberOfIds(); ++e){

		if (result -> GetId(e) != PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetQueriedPointId()){
			this -> neighbors_indices_vector.push_back(result -> GetId(e));
			
			double neighbor[3];
			this -> parent -> get_wrapped_shape_data()[name]-> get_polydata() -> GetPoint(result -> GetId(e),neighbor);

			double distance = std::sqrt(vtkMath::Distance2BetweenPoints(queried_point,neighbor));
			this -> neighbors_distances_vector.push_back(distance);


			vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
			sphere -> SetCenter (neighbor[0], neighbor[1],  neighbor[2]);

			sphere -> SetPhiResolution (20);
			sphere -> SetThetaResolution (20);
			sphere -> SetRadius (1);


			sphere -> Update();



			vtkSmartPointer<vtkDataSetMapper> neighborMapper = vtkSmartPointer<vtkDataSetMapper>::New();
			vtkSmartPointer<vtkActor> neighborActor = vtkSmartPointer<vtkActor>::New();
			neighborActor->SetMapper(neighborMapper);

			neighborMapper -> SetInputConnection(sphere -> GetOutputPort());
			neighborActor -> GetProperty() -> SetColor(1,1,0);
			neighborActor -> GetProperty() -> VertexVisibilityOff();
			neighborActor -> GetProperty() -> EdgeVisibilityOff();

			neighbors_actors_vector.push_back(neighborActor);

			this -> parent -> get_renderer() -> AddActor(neighborActor);

		}

	}

	this -> update_neighbors_positions();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}


void VertexEditionWindow::update_neighbors_positions(){

	double length = PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) ->  GetSelectedBodySize();
	double val_0 = 0.25 * length * (double(this -> position_slider -> value()) - 50)/50.;

	for (int e = 0; e < this -> neighbors_distances_vector.size(); ++e ){
		double n[3];

		PickInteractorStyle::SafeDownCast(this -> parent -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetNormalAtPoint(this -> neighbors_indices_vector[e],n);


		double val = std::exp(- std::pow(this -> neighbors_distances_vector[e]/(this -> neighbors_distances_vector.back()/3),2)) * val_0;
		
		double translation_vector[3] = { 
			val * n[0], 
			val * n[1],  
			val * n[2]
		};

		double position[3] = {
			translation_vector[0],
			translation_vector[1],
			translation_vector[2]
		};

		this -> neighbors_actors_vector[e] -> SetPosition(position);
		this -> neighbors_actors_vector[e] -> Modified();

	}

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

