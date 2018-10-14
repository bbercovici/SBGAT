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

#include <QScrollArea>

#include "ShapePropertiesWidget.hpp"
#include "SurfacePGMWindow.hpp"

#include <SBGATMassProperties.hpp>
#include <SBGATTrajectory.hpp>
#include <OrbitConversions.hpp>
#include <SBGATPolyhedronGravityModel.hpp>

using namespace SBGAT_GUI;

SurfacePGMWindow::SurfacePGMWindow(Mainwindow * parent) {

	this -> parent = parent;


	QVBoxLayout * window_layout = new QVBoxLayout(this);
	QWidget * select_shape_widget = new QWidget(this);
	QHBoxLayout * select_shape_layout = new QHBoxLayout(select_shape_widget);

	this -> compute_surface_pgm_button = new QPushButton("Compute Surface PGM",this);
	this -> primary_prop_combo_box = new QComboBox (this);
	this -> primary_shape_properties_widget = new ShapePropertiesWidget(this ,"Shape properties");

	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);


	select_shape_layout -> addWidget(new QLabel("Shape model",this));
	select_shape_layout -> addWidget(this -> primary_prop_combo_box);

	window_layout -> addWidget(select_shape_widget);
	window_layout -> addWidget(this -> primary_shape_properties_widget);
	window_layout -> addWidget(this -> compute_surface_pgm_button);
	window_layout -> addWidget(this -> button_box);



	this -> init();
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> compute_surface_pgm_button, SIGNAL(clicked()), this, SLOT(compute_surface_pgm()));


	window_layout -> addStretch(1);


}

void SurfacePGMWindow::init(){

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	if (wrapped_shape_data.size() == 0 ){
		this -> compute_surface_pgm_button -> setEnabled(false);
	}
	else{
		this -> compute_surface_pgm_button -> setEnabled(true);

		for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
			this -> primary_prop_combo_box -> insertItem(this -> primary_prop_combo_box -> count(),QString::fromStdString(it -> first));
		}
	}


}

void SurfacePGMWindow::compute_surface_pgm(){

	std::string selected_shape_name = this -> primary_prop_combo_box -> currentText().toStdString();
	vtkSmartPointer<vtkPolyData> selected_shape = this -> parent -> get_wrapped_shape_data()[selected_shape_name]-> get_polydata();

	vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();
	pgm_filter -> SetInputData(selected_shape);

	pgm_filter -> SetDensity(this -> primary_shape_properties_widget -> get_density());
	pgm_filter -> SetScaleMeters();
	pgm_filter -> Update();

	int numCells = selected_shape -> GetNumberOfCells();

	arma::vec::fixed<3> omega = this -> primary_shape_properties_widget -> get_spin();
	omega *= 2 * arma::datum::pi / this -> primary_shape_properties_widget -> get_period();


	std::vector<double> slopes,potentials,acc_magnitudes,acc_body_fixed_magnitudes;

    // The facets are browsed 
	for (int cellId=0; cellId < numCells; cellId++){

		vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
		ptIds -> Allocate(3);
		selected_shape -> GetCellPoints(cellId,ptIds);

		int p0_index = ptIds -> GetId(0);
		int p1_index = ptIds -> GetId(1);
		int p2_index = ptIds -> GetId(2);

		double p0[3];
		double p1[3];
		double p2[3];

		selected_shape -> GetPoint(p0_index,p0);
		selected_shape -> GetPoint(p1_index,p1);
		selected_shape -> GetPoint(p2_index,p2);

		arma::vec::fixed<3> p0_arma = {p0[0],p0[1],p0[2]};
		arma::vec::fixed<3> p1_arma = {p1[0],p1[1],p1[2]};
		arma::vec::fixed<3> p2_arma = {p2[0],p2[1],p2[2]};
		arma::vec::fixed<3> facet_center = 1./3 * (p0_arma + p1_arma + p2_arma);
		arma::vec::fixed<3> normal = arma::normalise(arma::cross(p1_arma - p0_arma,p2_arma - p0_arma));

		double potential,slope;
		arma::vec::fixed<3> acc,acc_body_fixed;
		
		pgm_filter -> GetPotentialAcceleration(facet_center,potential,acc);
		acc_body_fixed = acc - arma::cross(omega,arma::cross(omega,facet_center));
		slope = std::acos(arma::dot(-arma::normalise(acc_body_fixed),normal)) * 180./arma::datum::pi;

		slopes.push_back(slope);
		potentials.push_back(potential);
		acc_magnitudes.push_back(arma::norm(acc));
		acc_body_fixed_magnitudes.push_back(arma::norm(acc_body_fixed));

	}	

	std::shared_ptr<ModelDataWrapper> wrapper = this -> parent -> get_wrapped_shape_data()[selected_shape_name];

	wrapper -> set_potentials(potentials);
	wrapper -> set_acc_magnitudes(acc_magnitudes);
	wrapper -> set_slopes(slopes);
	wrapper -> set_acc_body_fixed_magnitudes(acc_body_fixed_magnitudes);
	wrapper -> get_mapper() -> ScalarVisibilityOff();

	this -> parent -> get_renderer() -> RemoveActor2D(wrapper -> get_colorbar_actor());

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

}



