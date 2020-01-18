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
#include "MassPropertiesWindow.hpp"
#include <SBGATMassProperties.hpp>
#include <SBGATMassPropertiesUQ.hpp>


using namespace SBGAT_GUI;

MassPropertiesWindow::MassPropertiesWindow(Mainwindow * parent) : AnalysesWindow(parent){

	this -> setWindowTitle("Compute Mass Properties");

	QWidget * button_widget = new QWidget(this);

	QHBoxLayout * button_widget_layout = new QHBoxLayout(button_widget);


	this -> compute_mass_properties_button = new QPushButton("Compute Mass Properties",this);

	this -> primary_shape_uncertainty_widget = new ShapeUncertaintyWidget(this ,"Shape uncertainty");

	// Creating the button box
	this -> button_box = new QDialogButtonBox(QDialogButtonBox::Ok);

	// Creating the output folder 
	this -> analyses_layout -> addWidget(this -> primary_shape_uncertainty_widget);

	this -> save_to_file_checkbox = new QCheckBox(this);

	QWidget * save_to_file_widget = new QWidget(this);
	QHBoxLayout * save_to_file_layout = new QHBoxLayout(save_to_file_widget);

	save_to_file_layout -> addWidget(new QLabel("Save To File",this));
	save_to_file_layout -> addWidget(this -> save_to_file_checkbox);

	this -> analyses_layout -> addWidget(save_to_file_widget);



	button_widget_layout -> addWidget(this -> compute_mass_properties_button);

	this -> analyses_layout -> addWidget(button_widget);

	this -> analyses_layout -> addWidget(this -> button_box);


	this -> init();
	connect(this -> button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(this -> compute_mass_properties_button, SIGNAL(clicked()), this, SLOT(compute_mass_properties()));
	connect(this -> prop_combo_box,SIGNAL(currentIndexChanged(int)),this -> primary_shape_uncertainty_widget,SLOT(clear()));



	this -> analyses_layout -> addStretch(1);


}

void MassPropertiesWindow::init(){

	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();	
	
	
	if (wrapped_shape_data.size() != 0 ){
		for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
			this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		}
	}

	

}

void MassPropertiesWindow::compute_mass_properties(){

	std::string selected_shape_name = this -> prop_combo_box -> currentText().toStdString();
	vtkSmartPointer<vtkPolyData> selected_shape = this -> parent -> get_wrapped_shape_data()[selected_shape_name]-> get_polydata();

	
	std::string opening_line = "### Computing surface PGM of " + selected_shape_name  + " ###";
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(opening_line));


	auto start = std::chrono::system_clock::now();

	vtkSmartPointer<SBGATMassProperties> mass_properties = vtkSmartPointer<SBGATMassProperties>::New();

	mass_properties -> SetInputData(selected_shape);
	mass_properties -> SetDensity(1.);// the computed mass properties do not depend on the density/mass, but only on the shape

	mass_properties -> SetScaleMeters();

	mass_properties -> Update();

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end-start;



	// An instance of SBGATMassPropertiesUQ is created to perform
	// uncertainty quantification from the PGM associated to the shape
	if (this -> primary_shape_uncertainty_widget -> get_active_tab_index() != 0){

		SBGATMassPropertiesUQ mass_uq;
		mass_uq.SetModel(mass_properties);
		mass_uq.PrecomputeMassPropertiesPartials();

		// The shape covariance is extracted
		arma::mat P_CC;

		QString cov_path;

		switch(this -> primary_shape_uncertainty_widget -> get_active_tab_index()){
			case 1:
			// The shape covariance is loaded from a file
			// No regularization is applied
			P_CC.load(this -> primary_shape_uncertainty_widget -> get_covariance_input_file());
			if (P_CC.n_rows != 3 * mass_properties -> GetN_vertices() || P_CC.n_cols != 3 * mass_properties -> GetN_vertices()){
				QMessageBox::warning(this, "Compute Mass Properties","The provided covariance matrix's dimensions do not match that of the considered shape. No output was produced");
				return;
			}
			else{
				mass_uq.SetCovariance(P_CC);
			}
			break;
			case 2:
			// The shape covariance is created from a global uncertainty description
			
			if (this -> primary_shape_uncertainty_widget -> get_save_global_covariance_to_file_checkbox()){
				cov_path = QFileDialog::getSaveFileName(this, tr("Save Shape Covariance To File"),"~/",tr("Text file (*.txt)"));
			}

			mass_uq.ComputeVerticesCovarianceGlobal(this -> primary_shape_uncertainty_widget -> get_global_sigma(),
				this -> primary_shape_uncertainty_widget -> get_global_correlation_distance());
			

			for (unsigned int k = 0; k < this -> primary_shape_uncertainty_widget -> get_global_covariance_regularization_number(); ++k){
				int regularized_eigenvalues = mass_uq.RegularizeCovariance();
				this -> parent -> log_console -> appendPlainText(QString::fromStdString("- " + std::to_string(regularized_eigenvalues) + " P_CC eigenvalues were regularized"));
				if (regularized_eigenvalues == 0) 
					break;
			}


			break;
			case 3:
			// The shape covariance is created from a local uncertainty description

			if (this -> primary_shape_uncertainty_widget -> local_uncertainty_table_widget -> rowCount() == 0){
				QMessageBox::warning(this, "Compute Mass Properties","There are no uncertainty regions for this shape");
				this -> parent -> log_console -> appendPlainText(QString::fromStdString("The surface PGM computation was interrupted"));
				
				return;
				break;
			}
			else{
				if (this -> primary_shape_uncertainty_widget -> get_save_local_covariance_to_file_checkbox()){
					cov_path = QFileDialog::getSaveFileName(this, tr("Save Shape Covariance To File"),"~/",tr("Text file (*.txt)"));
				}
				for (int i = 0; i < this -> primary_shape_uncertainty_widget -> local_uncertainty_table_widget -> rowCount(); ++i){
					

					int region_center = qobject_cast<QSpinBox*>(this -> primary_shape_uncertainty_widget -> local_uncertainty_table_widget -> cellWidget(i,0)) -> value();
					double error_standard_dev = qobject_cast<QDoubleSpinBox*>(this -> primary_shape_uncertainty_widget -> local_uncertainty_table_widget -> cellWidget(i,1)) -> value();
					double correlation_distance = qobject_cast<QDoubleSpinBox*>(this -> primary_shape_uncertainty_widget -> local_uncertainty_table_widget -> cellWidget(i,2)) -> value();

					mass_uq.AddRadialUncertaintyRegionToCovariance(region_center,error_standard_dev,correlation_distance);

				}

				for (unsigned int k = 0; k < this -> primary_shape_uncertainty_widget -> get_local_covariance_regularization_number(); ++k){
					int regularized_eigenvalues = mass_uq.RegularizeCovariance();
					this -> parent -> log_console -> appendPlainText(QString::fromStdString("- " + std::to_string(regularized_eigenvalues) + " P_CC eigenvalues were regularized"));
					if (regularized_eigenvalues == 0) 
						break;
				}
			}

			
			
		}

		// 

		double volume_variance;
		arma::mat::fixed<3,3> cov_cm;
		arma::mat::fixed<6,6> cov_I;
		arma::mat::fixed<3,3> cov_sigma;

		mass_uq.GetMassPropertiesUncertainties(volume_variance,cov_cm,cov_I,cov_sigma);


		if (cov_path.length() > 0){
			mass_uq.GetVerticesCovariance().save(cov_path.toStdString(),arma::raw_ascii);
		}

	}
	else if (this -> save_to_file_checkbox -> isChecked()){

		std::string default_name;

		std::string name = this -> prop_combo_box -> currentText().toStdString();
		auto shape_data = this -> parent -> get_wrapped_shape_data();

		if ( shape_data.find(name)!= shape_data.end()){
			default_name = "./" + name + "_mass_properties.json";
		}

		QString path = QFileDialog::getSaveFileName(this, tr("Save Mass Properties To File"),
			QString::fromStdString(default_name),
			tr("JSON file (*.json)"));
		if (path.size() > 0){

			mass_properties -> SaveMassProperties(path.toStdString());
		}
	}

	std::stringstream ss;
	ss.str(std::string());
	ss.precision(10);

	opening_line = "### Computing geometric measures of " + selected_shape_name + " ###";
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(opening_line));


	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Surface area (m^2) :"));
	this -> parent -> log_console -> appendPlainText(" " + QString::number(mass_properties -> GetSurfaceArea ()));

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Volume (m^3) :"));
	this -> parent -> log_console -> appendPlainText(" " + QString::number(mass_properties -> GetVolume()));


	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Average radius (m) :"));
	this -> parent -> log_console -> appendPlainText(" " + QString::number(mass_properties -> GetAverageRadius()));

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Bounding box (m) :"));

	double xmin,xmax, ymin, ymax, zmin, zmax;
	mass_properties -> GetBoundingBox(xmin,xmax,ymin,ymax,zmin,zmax);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("-- Min: " + std::to_string(xmin) + " "+ std::to_string(ymin) + " "+ std::to_string(zmin)));
	this -> parent -> log_console -> appendPlainText(QString::fromStdString("-- Max: " + std::to_string(xmax) + " "+ std::to_string(ymax) + " "+ std::to_string(zmax)));

	ss.str(std::string());
	ss.precision(10);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Center of mass (m) :"));
	mass_properties -> GetCenterOfMass().t().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

	ss.str(std::string());
	ss.precision(10);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Normalized inertia tensor :"));
	mass_properties -> GetNormalizedInertiaTensor().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

	ss.str(std::string());
	ss.precision(10);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Normalized inertia moments :"));
	mass_properties -> GetNormalizedInertiaMoments().t().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

	ss.str(std::string());
	ss.precision(10);


	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Unit-density inertia tensor (m^5) :"));
	mass_properties -> GetUnitDensityInertiaTensor().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));


	ss.str(std::string());
	ss.precision(10);


	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Unit-density inertia moments (m^5) :"));
	mass_properties -> GetUnitDensityInertiaMoments().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

	ss.str(std::string());
	ss.precision(10);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\b- Principal dimensions (m) :"));
	mass_properties -> GetPrincipalDimensions().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

	ss.str(std::string());
	ss.precision(10);

	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\b- Principal axes :"));
	mass_properties -> GetPrincipalAxes().raw_print(ss);
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(ss.str()));


	this -> parent -> log_console -> appendPlainText(QString::fromStdString("\n- Done computing in ")
		+ QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

	std::string closing_line(opening_line.length() - 1, '#');
	closing_line.append("\n");
	this -> parent -> log_console -> appendPlainText(QString::fromStdString(closing_line));











}


