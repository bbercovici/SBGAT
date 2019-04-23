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
#include "ShapeUncertaintyWidget.hpp"
#include "SurfacePGMWindow.hpp"

#include <SBGATMassProperties.hpp>
#include <SBGATTrajectory.hpp>
#include <OrbitConversions.hpp>
#include <SBGATPolyhedronGravityModel.hpp>
#include <SBGATPolyhedronGravityModelUQ.hpp>
using namespace SBGAT_GUI;

SurfacePGMWindow::SurfacePGMWindow(Mainwindow * parent) {

	this -> parent = parent;
	this -> setWindowTitle("Compute/Load Surface PGM");

	QVBoxLayout * window_layout = new QVBoxLayout(this);
	QWidget * select_shape_widget = new QWidget(this);
	QWidget * button_widget = new QWidget(this);

	QHBoxLayout * select_shape_layout = new QHBoxLayout(select_shape_widget);
	QHBoxLayout * button_widget_layout = new QHBoxLayout(button_widget);


	this -> compute_surface_pgm_button = new QPushButton("Compute Surface PGM",this);
	this -> load_surface_pgm_button = new QPushButton("Load Surface PGM",this);

	this -> primary_prop_combo_box = new QComboBox (this);
	this -> primary_shape_properties_widget = new ShapePropertiesWidget(this ,"Shape properties");
	this -> primary_shape_uncertainty_widget = new ShapeUncertaintyWidget(this ,"Shape uncertainty");


	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	// Creating the output folder 
	this ->  open_output_file_dialog_button = new QPushButton("Select output file",this);


	select_shape_layout -> addWidget(new QLabel("Shape model",this));
	select_shape_layout -> addWidget(this -> primary_prop_combo_box);

	window_layout -> addWidget(select_shape_widget);

	window_layout -> addWidget(this -> primary_shape_properties_widget);
	window_layout -> addWidget(this -> primary_shape_uncertainty_widget);

	
	button_widget_layout -> addWidget(this -> compute_surface_pgm_button);
	button_widget_layout -> addWidget(this -> load_surface_pgm_button);

	window_layout -> addWidget(this -> open_output_file_dialog_button);
	window_layout -> addWidget(button_widget);

	window_layout -> addWidget(this -> button_box);



	this -> compute_surface_pgm_button -> setEnabled(false);

	this -> init();
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> compute_surface_pgm_button, SIGNAL(clicked()), this, SLOT(compute_surface_pgm()));
	connect(this -> load_surface_pgm_button, SIGNAL(clicked()), this, SLOT(load_surface_pgm()));
	connect(this -> open_output_file_dialog_button,SIGNAL(clicked()),this,
		SLOT(open_output_file_dialog()));


	window_layout -> addStretch(1);


}

void SurfacePGMWindow::init(){

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	
	if (wrapped_shape_data.size() != 0 ){
		for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
			this -> primary_prop_combo_box -> insertItem(this -> primary_prop_combo_box -> count(),QString::fromStdString(it -> first));
		}
	}

	if(wrapped_shape_data.size() == 0){
		this -> load_surface_pgm_button -> setEnabled(false);
	}

}

void SurfacePGMWindow::compute_surface_pgm(){

	std::string selected_shape_name = this -> primary_prop_combo_box -> currentText().toStdString();
	vtkSmartPointer<vtkPolyData> selected_shape = this -> parent -> get_wrapped_shape_data()[selected_shape_name]-> get_polydata();

	arma::vec::fixed<3> omega = this -> primary_shape_properties_widget -> get_spin();

	if (this -> primary_shape_properties_widget -> get_period() == 0){
		QMessageBox::warning(this, "Evaluate Surface PGM", "The rotation period must be strictly greater than 0!");
		return;
	}

	std::string opening_line = "### Computing surface PGM of " + selected_shape_name  + " ###";
	this -> parent ->log_console -> appendPlainText(QString::fromStdString(opening_line));

	omega *= 2 * arma::datum::pi / this -> primary_shape_properties_widget -> get_period();

	std::vector<double> slopes,
	inertial_potentials,
	body_fixed_potentials,
	inertial_acc_magnitudes,
	body_fixed_acc_magnitudes,
	analytical_slope_sds;



	int numCells = selected_shape -> GetNumberOfCells();
	std::vector<unsigned int> queried_elements;
	
	for (unsigned int i = 0; i < static_cast<unsigned int>(numCells); ++i){
		queried_elements.push_back(i);
	}

	auto start = std::chrono::system_clock::now();


	SBGATPolyhedronGravityModel::ComputeSurfacePGM(selected_shape,
		queried_elements,
		true,
		this -> primary_shape_properties_widget -> get_density(),
		omega,
		slopes,
		inertial_potentials,
		body_fixed_potentials,
		inertial_acc_magnitudes,
		body_fixed_acc_magnitudes);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;
	


	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();

	mass_properties -> SetInputData(selected_shape);
	mass_properties -> Update();
	double mass = mass_properties -> GetVolume() * this -> primary_shape_properties_widget -> get_density();







	


	// An instance of SBGATPolyhedronGravityModelUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	if (this -> primary_shape_uncertainty_widget -> get_active_tab_index() != 0){


		// An instance of SBGATPolyhedronGravityModel is created to evaluate the PGM of 
	// the considered polytdata
		vtkSmartPointer<SBGATPolyhedronGravityModel> pgm_filter = vtkSmartPointer<SBGATPolyhedronGravityModel>::New();

		pgm_filter -> SetInputData(selected_shape);
		pgm_filter -> SetDensity(this -> primary_shape_properties_widget -> get_density());

		pgm_filter -> SetScaleMeters();
		pgm_filter -> SetOmega(omega);
		pgm_filter -> Update();


		SBGATPolyhedronGravityModelUQ pgm_uq;
		pgm_uq.SetModel(pgm_filter);
		pgm_uq.SetPeriodErrorStandardDeviation(0);
		pgm_uq.PrecomputeMassPropertiesPartials();
		std::vector<double> analytical_slope_variances;


		// The shape covariance is extracted
		bool early_exit = false;
		arma::mat P_CC;

		switch(this -> primary_shape_uncertainty_widget -> get_active_tab_index()){
			case 1:
			// The shape covariance is loaded from a file
			P_CC.load(this -> primary_shape_uncertainty_widget -> get_covariance_input_file());
			if (P_CC.n_rows != 3 * pgm_filter -> GetN_vertices() || P_CC.n_cols != 3 * pgm_filter -> GetN_vertices()){
				QMessageBox::warning(this, "Compute Surface PGM","The provided covariance matrix's dimensions do not match that of the considered shape. No output was produced");
				early_exit = true;
			}
			else{
				pgm_uq.SetCovariance(P_CC);
			}
			break;
			case 2:
			// The shape covariance is created from a global uncertainty description
			pgm_uq.ComputeVerticesCovarianceGlobal(this -> primary_shape_uncertainty_widget -> get_global_sigma(),
				this -> primary_shape_uncertainty_widget -> get_global_correlation_distance());
			

			for (unsigned int k = 0; k < this -> primary_shape_uncertainty_widget -> get_global_covariance_regularization_number(); ++k){
				int regularized_eigenvalues = pgm_uq.RegularizeCovariance();
				this -> parent -> log_console -> appendPlainText(QString::fromStdString("- " + std::to_string(regularized_eigenvalues) + " P_CC eigenvalues were regularized"));
				if (regularized_eigenvalues == 0) 
					break;
			}


			break;
			case 3:
			QMessageBox::warning(this, "Compute Surface PGM","Not implemented yet");
			early_exit = true;
			break;
		}

		if (early_exit){
			return;
		}


		pgm_uq.GetVarianceSlopes(analytical_slope_variances,queried_elements,false);



		for (auto it = analytical_slope_variances.begin(); it != analytical_slope_variances.end(); ++it){
			analytical_slope_sds.push_back(std::sqrt(*it));
		}

		SBGATPolyhedronGravityModel::SaveSurfacePGM(selected_shape,
			queried_elements,
			true,
			mass,
			omega,
			slopes,
			inertial_potentials,
			body_fixed_potentials,
			inertial_acc_magnitudes,
			body_fixed_acc_magnitudes,
			analytical_slope_sds,
			this -> output_path );




		QString cov_path = QFileDialog::getSaveFileName(this, tr("Save Shape Covariance To File"),
			"~/",
			tr("Text file (*.txt)"));

		if (cov_path.length() > 0){
			pgm_uq.GetVerticesCovariance().save(cov_path.toStdString(),arma::raw_ascii);
		}


	}
	else{

		SBGATPolyhedronGravityModel::SaveSurfacePGM(selected_shape,
			queried_elements,
			true,
			mass,
			omega,
			slopes,
			inertial_potentials,
			body_fixed_potentials,
			inertial_acc_magnitudes,
			body_fixed_acc_magnitudes,
			this -> output_path );
	}


	std::shared_ptr<ModelDataWrapper> wrapper = this -> parent -> get_wrapped_shape_data()[selected_shape_name];

	wrapper -> set_inertial_potentials(inertial_potentials);
	wrapper -> set_body_fixed_potentials(body_fixed_potentials);

	wrapper -> set_inertial_acc_magnitudes(inertial_acc_magnitudes);
	wrapper -> set_body_fixed_acc_magnitudes(body_fixed_acc_magnitudes);
	wrapper -> set_slopes(slopes);

	if (this -> primary_shape_uncertainty_widget -> get_active_tab_index() != 0){
		wrapper -> set_slope_sds(analytical_slope_sds);
	}

	wrapper -> get_mapper() -> ScalarVisibilityOff();

	this -> parent -> get_renderer() -> RemoveActor2D(wrapper -> get_colorbar_actor());

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

	std::string displayed_line = "- Saved surface-evaluated PGM of " + selected_shape_name + " to " + this -> output_path;

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("- Done computing surface PGM in " + std::to_string(elapsed_seconds.count()) +  " seconds."));
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(displayed_line));

	std::string closing_line(opening_line.length() - 1, '#');

	closing_line.append("\n");

	this -> parent -> log_console -> appendPlainText(QString::fromStdString(closing_line));










}


void SurfacePGMWindow::open_output_file_dialog(){

	std::string default_name;

	std::string name = this -> primary_prop_combo_box -> currentText().toStdString();
	auto shape_data = this -> parent -> get_wrapped_shape_data();
	
	if ( shape_data.find(name)!= shape_data.end()){
		default_name = "./" + name + "_surface_pgm.json";
	}

	QString path = QFileDialog::getSaveFileName(this, tr("Save Surface PGM To File"),
		QString::fromStdString(default_name),
		tr("JSON file (*.json)"));

	this -> output_path = path.toStdString();

	if (this -> output_path.size() > 0){
		this -> compute_surface_pgm_button -> setEnabled(true);

	}
}

void SurfacePGMWindow::load_surface_pgm(){
	std::string selected_shape_name = this -> primary_prop_combo_box -> currentText().toStdString();

	QString path = QFileDialog::getOpenFileName(this, tr("Open Surface PGM File"),
		"~",
		tr("JSON (*.json)"));
	
	if (path.size() == 0){
		return;
	}

	std::vector<double> slopes,inertial_potentials,
	body_fixed_potentials,
	inertial_acc_magnitudes,
	body_fixed_acc_magnitudes,
	slope_sds;
	double mass;
	arma::vec::fixed<3> omega;

	try{
		std::shared_ptr<ModelDataWrapper> wrapper = this -> parent -> get_wrapped_shape_data()[selected_shape_name];

		try{
			SBGATPolyhedronGravityModel::LoadSurfacePGM(mass,
				omega,
				slopes,
				inertial_potentials,
				body_fixed_potentials,
				inertial_acc_magnitudes,
				body_fixed_acc_magnitudes,
				slope_sds,
				path.toStdString());
		}

		catch(std::runtime_error & e){
			QMessageBox::warning(this, "Load Surface PGM", e.what());
			return;
		}

		if (slopes.size() != wrapper -> get_polydata() -> GetNumberOfCells()){
			QMessageBox::warning(this, 
				"Load Surface PGM", QString::fromStdString("Error: the loaded surface PGM (" + std::to_string(slopes.size())+ 
					" facets) does not match the selected shape resolution (" + std::to_string(wrapper -> get_polydata() -> GetNumberOfCells()) + " facets)"));
			return;
		}

		wrapper -> set_slopes(slopes);
		wrapper -> set_slope_sds(slope_sds);


		wrapper -> set_inertial_potentials(inertial_potentials);
		wrapper -> set_body_fixed_potentials(body_fixed_potentials);

		wrapper -> set_inertial_acc_magnitudes(inertial_acc_magnitudes);
		wrapper -> set_body_fixed_acc_magnitudes(body_fixed_acc_magnitudes);
		
		wrapper -> get_mapper() -> ScalarVisibilityOff();

		this -> parent -> get_renderer() -> RemoveActor2D(wrapper -> get_colorbar_actor());

		this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
		this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Done loading surface PGM from file " + path.toStdString()));

	}
	catch(std::runtime_error & e ){
		QMessageBox::warning(this, "Load Surface PGM", "An error occured loading the Surface PGM file.");
		return;
	}

}


