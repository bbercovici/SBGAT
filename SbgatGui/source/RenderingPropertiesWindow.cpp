#include "RenderingPropertiesWindow.hpp"


using namespace SBGAT_GUI;

RenderingPropertiesWindow::RenderingPropertiesWindow(Mainwindow * parent) : QDialog(parent,Qt::WindowStaysOnTopHint) {

	this -> parent = parent;
	this -> setWindowTitle("Rendering properties");

	QVBoxLayout * main_layout = new QVBoxLayout(this);
	QGroupBox * focus_prop_group = new QGroupBox(tr("Camera focus"));
	QGroupBox * light_group = new QGroupBox(tr("Lights"));
	QGroupBox * current_light_group = new QGroupBox(tr("Current light"));
	QGroupBox * new_light_group = new QGroupBox(tr("New light"));
	QGroupBox * shadow_group = new QGroupBox(tr("Shadows"));

	QPushButton * add_light_button = new QPushButton(tr("Add light"),this);
	QCheckBox * enable_mutual_shadows_check = new QCheckBox(this);

	
	this -> remove_light_button = new QPushButton(tr("Remove light"),this);
	this -> new_light_combo_box = new QComboBox(this);
	this -> current_light_combo_box = new QComboBox(this);
	this -> prop_combo_box = new QComboBox (this);

	QLabel * focus_prop_label = new QLabel("Select prop to focus on",this);
	QLabel * shadow_label = new QLabel("Enable mutual shadows",this);


	QGridLayout * focus_prop_layout = new QGridLayout(focus_prop_group);
	QGridLayout * shadow_group_layout = new QGridLayout(shadow_group);


	focus_prop_layout -> addWidget(focus_prop_label, 0, 0, 1, 1);
	focus_prop_layout -> addWidget(prop_combo_box, 0, 1, 1, 1);


	QVBoxLayout * light_group_layout = new QVBoxLayout(light_group);
	QGridLayout * current_light_group_layout = new QGridLayout(current_light_group);
	QGridLayout * new_light_group_layout = new QGridLayout(new_light_group);


	light_group_layout -> addWidget(current_light_group);
	light_group_layout -> addWidget(new_light_group);	


	current_light_group_layout -> addWidget(this -> current_light_combo_box, 0, 0, 1, 1);
	current_light_group_layout -> addWidget(this -> remove_light_button, 0, 1, 1, 1);
	new_light_group_layout -> addWidget(this -> new_light_combo_box, 0, 0 , 1, 1);
	new_light_group_layout -> addWidget(add_light_button, 0, 1, 1, 1);

	focus_prop_layout -> addWidget(focus_prop_label, 0, 0, 1, 1);
	focus_prop_layout -> addWidget(prop_combo_box, 0, 1, 1, 1);

	shadow_group_layout -> addWidget(shadow_label,0,0,1,1);
	shadow_group_layout -> addWidget(enable_mutual_shadows_check,0,1,1,1);

	main_layout -> addWidget(focus_prop_group);
	main_layout -> addWidget(light_group);
	main_layout -> addWidget(shadow_group);


	QDialogButtonBox * button_box = new QDialogButtonBox(QDialogButtonBox::Ok
		| QDialogButtonBox::Cancel);


	connect(add_light_button,SIGNAL(clicked()),this,SLOT(add_light()));
	connect(this -> remove_light_button,SIGNAL(clicked()),this,SLOT(remove_light()));


	connect(button_box, SIGNAL(accepted()), this, SLOT(accept()));
	connect(button_box, SIGNAL(rejected()), this, SLOT(close()));
	connect(this -> prop_combo_box,SIGNAL(currentIndexChanged(int)),this,SLOT(change_focus()));
	connect(enable_mutual_shadows_check,SIGNAL(stateChanged(int)), this , SLOT(enable_mutual_shadows(int)));
	
	this -> init();

}

void RenderingPropertiesWindow::init(){
	
	std::vector<std::string> light_type = {"Scene light","Headlight","Camera light"};


	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();
	
	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		
	}

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),QString::fromStdString(it -> first));
		
	}


	auto lights = this -> parent -> get_renderer() -> GetLights();
	lights -> InitTraversal();

	for(vtkIdType i = 0; i < lights -> GetNumberOfItems(); i++){
		auto light = lights -> GetNextItem();
		this -> current_light_combo_box -> insertItem(this -> current_light_combo_box -> count(),QString::fromStdString(light_type[light -> GetLightType()]));
	}


	this -> new_light_combo_box -> insertItem(0,tr("Scene light"));
	this -> new_light_combo_box -> insertItem(1,tr("Headlight"));
	this -> new_light_combo_box -> insertItem(2,tr("Camera light"));

	this -> make_light_box_consistent();

}


void RenderingPropertiesWindow::prop_removed_slot(){

	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

	for (int i = 0; i < this -> prop_combo_box -> count(); ++i){
		QString prop_name = this -> prop_combo_box -> itemText(i);
		if (wrapped_spacecraft_data.find(prop_name.toStdString()) == wrapped_spacecraft_data.end()&& wrapped_shape_data.find(prop_name.toStdString()) == wrapped_shape_data.end()){
			this -> prop_combo_box -> removeItem(i);
			break;

		}

	}




}







void RenderingPropertiesWindow::prop_added_slot(){
	auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
	auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

	for (auto it = wrapped_spacecraft_data.begin(); it != wrapped_spacecraft_data.end(); ++it){
		if (this -> prop_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}

	for (auto it = wrapped_shape_data.begin(); it != wrapped_shape_data.end(); ++it){
		if (this -> prop_combo_box  -> findText(QString::fromStdString(it -> first)) == -1){
			this -> prop_combo_box -> insertItem(this -> prop_combo_box -> count(),
				QString::fromStdString(it -> first));
			break;
		}
	}



}

void RenderingPropertiesWindow::change_focus(){

	if (this -> prop_combo_box -> count() > 0 ){

		auto wrapped_spacecraft_data = this -> parent -> get_wrapped_spacecraft_data();
		auto wrapped_shape_data = this -> parent -> get_wrapped_shape_data();

		std::string current_prop_name = this -> prop_combo_box -> currentText().toStdString();

		vtkSmartPointer<vtkActor> prop_to_focus_on;

		if (wrapped_spacecraft_data.find(current_prop_name) != wrapped_spacecraft_data.end()){
			prop_to_focus_on = wrapped_spacecraft_data[current_prop_name] -> get_actor();
		}
		else{
			prop_to_focus_on = wrapped_shape_data[current_prop_name] -> get_actor();
		}



		vtkSmartPointer<vtkCamera> camera = this -> parent -> get_renderer() -> GetActiveCamera();
		camera -> SetFocalPoint(prop_to_focus_on -> GetPosition());

		this -> parent -> get_renderer() -> Modified();
		this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();
	}


}

void RenderingPropertiesWindow::add_light(){

	std::vector<std::string> light_type = {"Scene light","Headlight","Camera light"};

	vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
	vtkSmartPointer<vtkLightCollection> lights = this -> parent -> get_renderer() -> GetLights();


	light -> SetFocalPoint(this -> parent ->get_renderer() -> GetActiveCamera() -> GetFocalPoint());
	light -> SetPosition(this -> parent ->get_renderer() -> GetActiveCamera() -> GetPosition());
	light -> SetColor(1,1,1);
	light -> SetIntensity(1);
	light -> SetLightType(this -> new_light_combo_box -> currentIndex()); // 0: scene , 1: headlight , 2: camera light

	this -> parent -> get_renderer() -> AddLight(light);

	
	auto N_lights = this -> current_light_combo_box -> count();
	this -> current_light_combo_box -> insertItem(N_lights,QString::fromStdString(light_type[light -> GetLightType()]));
	this -> current_light_combo_box -> setCurrentIndex(N_lights);
	this -> make_light_box_consistent();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}

void RenderingPropertiesWindow::remove_light(){
	int current_light = this -> current_light_combo_box -> currentIndex();
	vtkSmartPointer<vtkLightCollection> lights = this -> parent -> get_renderer() -> GetLights();
	lights -> InitTraversal();

	
	for(vtkIdType i = 0; i < lights -> GetNumberOfItems(); i++){
		vtkLight * light = lights -> GetNextItem();
		if (i == current_light ){
			this -> parent -> get_renderer() -> RemoveLight(light);
			this -> current_light_combo_box -> removeItem(current_light);
			
			this -> current_light_combo_box -> setCurrentIndex(this -> current_light_combo_box -> count() - 1);
			

			break;
		}
	}

	this -> make_light_box_consistent();

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();

	
}


void RenderingPropertiesWindow::make_light_box_consistent(){
	auto lights = this -> parent -> get_renderer() -> GetLights();

	if (lights -> GetNumberOfItems() == 1) {
		this -> remove_light_button -> setDisabled(true);
	}

	else{
		this -> remove_light_button -> setEnabled(true);
	}
	this -> current_light_combo_box -> repaint();
	this -> remove_light_button -> repaint();




}

void RenderingPropertiesWindow::enable_mutual_shadows(int state){

	vtkSmartPointer<vtkShadowMapPass> shadows = vtkSmartPointer<vtkShadowMapPass>::New();

	vtkSmartPointer<vtkSequencePass> seq = vtkSmartPointer<vtkSequencePass>::New();
	vtkSmartPointer<vtkRenderPassCollection> passes = vtkSmartPointer<vtkRenderPassCollection>::New();
	
	passes -> AddItem(shadows -> GetShadowMapBakerPass());
	passes -> AddItem(shadows);
	seq -> SetPasses(passes);

	vtkSmartPointer<vtkCameraPass> cameraP = vtkSmartPointer<vtkCameraPass>::New();
	cameraP->SetDelegatePass(seq);

	vtkOpenGLRenderer *glrenderer = vtkOpenGLRenderer::SafeDownCast(this -> parent -> get_renderer());
	glrenderer->SetPass(cameraP);

	this -> parent -> qvtkWidget -> GetRenderWindow() -> Render();


}
