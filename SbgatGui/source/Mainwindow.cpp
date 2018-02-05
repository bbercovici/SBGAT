#include "Mainwindow.hpp"
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFileDialog>
#include <QInputDialog>
#include <QTextStream>
#include <QMessageBox>
#include <QRegularExpression>
#include <QStringList>
#include <QPushButton>
#include <QThread>
#include <QHeaderView>

#include <vtkAreaPicker.h>
#include <vtkAxesActor.h>
#include <vtkPolygon.h>
#include <vtkPolyDataMapper.h>
#include <vtkDoubleArray.h>
#include <vtkScalarBarActor.h>
#include <vtkLookupTable.h>
#include <vtkCamera.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkParametricSpline.h>
#include <vtkLightActor.h>
#include <vtkLight.h>
#include <vtkLightCollection.h>
#include <vtkParametricFunctionSource.h>
#include <vtkOBJReader.h>
#include <vtkCenterOfMass.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkActor2DCollection.h>
#include <vtkCellData.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <ShapeModel.hpp>
#include <DynamicAnalyses.hpp>


#include "SettingsWindow.hpp"
#include "MoveAlongTrajectoryWindow.hpp"
#include "RenderingPropertiesWindow.hpp"
#include "SBGATMassProperties.hpp"
#include "Worker.hpp"


using namespace SBGAT_GUI;

// shortcut to interactor modes
#define INTERACTOR_IS_ORIENT 0
#define INTERACTOR_IS_SELECT 1


Mainwindow::Mainwindow() {

    // The GUI elements are created
    this -> setupUi();

    // The elements interfacing with SbgatCore are created
    this -> frame_graph = std::make_shared<SBGAT_CORE::FrameGraph>();
    this -> frame_graph -> add_frame("inertial_default"); // a default inertial frame of reference is created

}

void Mainwindow::setupUi() {
    this -> resize(1024, 768);

    // The widget are created
    this -> left_dockwidget = new QDockWidget(this);
    this -> right_dockwidget = new QDockWidget(this);
    this -> qvtkWidget = new QVTKOpenGLWidget(this);
    this -> status_bar = new QStatusBar(this);
    this -> log_console = new QPlainTextEdit(this);
    this -> log_console -> setReadOnly(true);
    this -> prop_table = new QTableWidget(0, 4, this);


    // The status bar is populated
    this -> setStatusBar(this -> status_bar);
    this -> statusBar() -> showMessage("Ready");

    // Headers are added to the shape table
    QStringList header_lists = {"Name", "Show", "Properties" , "Erase"};
    this -> prop_table -> setHorizontalHeaderLabels(header_lists);
    this -> prop_table -> horizontalHeader()->setStretchLastSection(true);
    this -> prop_table -> setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    
    // Selecting an item in the table highlights the entire row
    this -> prop_table -> setSelectionBehavior(QAbstractItemView::SelectRows);
    this -> prop_table -> setSelectionMode(QAbstractItemView::SingleSelection);


    // Central window
    this -> setCentralWidget(qvtkWidget);
    this -> setWindowTitle(QStringLiteral("SBGAT (WIP)"));

    // Actions and menus are created
    this -> createActions();
    this -> createMenus();
    this -> update_actions_availability();

    // The rendering window and its props are initialized
    this -> init_rendering_window();

    // A slot is connected to the signal sent by the table when a new selection is made
    connect(this -> prop_table, SIGNAL(currentItemChanged(QTableWidgetItem * , QTableWidgetItem * )),
        this, SLOT(update_GUI_changed_prop()));

    this -> qvtkWidget -> update();


    this -> qvtkWidget -> GetRenderWindow() -> Render();

    // The lateral dockwidgets are initialized
    // This is delayed until after the renderer is updated so
    // that the default props (default light...) have been instantiated

    this -> show();

    this -> init_left_dockwidget();
    this -> init_right_dockwidget();



}

void Mainwindow::init_rendering_window(){
 // A VTK renderer is created and linked with the qvtk widget
    this -> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> render_window = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
    this -> qvtkWidget -> SetRenderWindow(render_window);
    this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(this -> renderer);


    this -> renderer -> SetGradientBackground (true);
    this -> renderer -> SetBackground (0.5, 0.5, 1);

    vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();

    vtkSmartPointer<vtkInteractorStyleSwitch> style =
    vtkSmartPointer<vtkInteractorStyleSwitch>::New();

    render_window -> GetInteractor() -> SetInteractorStyle( style );
    render_window -> GetInteractor() -> SetPicker(areaPicker);

    vtkSmartPointer<vtkAxesActor> axes =
    vtkSmartPointer<vtkAxesActor>::New();

    this -> orientation_widget =
    vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    orientation_widget -> SetOrientationMarker( axes );
    orientation_widget -> SetInteractor(render_window -> GetInteractor() );
    orientation_widget -> SetViewport( 0.0, 0.0, 0.2, 0.2 );
    orientation_widget -> SetEnabled( 1 );
    orientation_widget -> InteractiveOff();

}



void Mainwindow::init_right_dockwidget(){

 // The lateral dockwidget is filled up
    this -> right_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );

    QWidget * right_dockwidget_container = new QWidget(this);
    QVBoxLayout * right_dockwidget_container_layout = new QVBoxLayout();

    right_dockwidget_container -> setLayout(right_dockwidget_container_layout);
    right_dockwidget_container_layout -> addWidget( this -> prop_table );
    right_dockwidget_container_layout -> addWidget( this -> log_console );

    this -> right_dockwidget -> setWidget(right_dockwidget_container);
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> right_dockwidget);
    this -> right_dockwidget -> hide();

    // VTK adds a Head light by default to the renderer
    std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();
    vtkSmartPointer<vtkLightCollection> lights = this -> get_renderer() -> GetLights();
    vtkSmartPointer<vtkLightActor> lightActor = vtkSmartPointer<vtkLightActor>::New();
    
    lights -> InitTraversal();

    vtkLight * light = lights -> GetNextItem();

    lightActor -> SetLight(light);    
    lightActor -> SetVisibility(1);
    model_data -> set_light(light);
    model_data -> set_light_actor(lightActor);
    std::string name = "Head light";

    this -> wrapped_light_data[name] = model_data;
    this -> get_renderer() -> AddViewProp(lightActor);

    // The prop table is updated to show the default light
    this -> add_prop_to_table_widget(name);


}


void Mainwindow::init_left_dockwidget(){

     // The lateral dockwidget is filled up
    this -> left_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );

    QWidget * left_dockwidget_container = new QWidget(this);
    QVBoxLayout * left_dockwidget_container_layout = new QVBoxLayout();

    QGroupBox * light_creation_group = new QGroupBox(tr("Add Rendering Lights"));
    QVBoxLayout * light_creation_layout = new QVBoxLayout(light_creation_group);
    QPushButton * add_scene_light_button = new QPushButton(tr("Scene light"),this);
    QPushButton * add_head_light_button = new QPushButton(tr("Head light"),this);
    QPushButton * add_camera_light_button = new QPushButton(tr("Camera light"),this);

    light_creation_layout -> addWidget(add_scene_light_button);
    light_creation_layout -> addWidget(add_head_light_button);
    light_creation_layout -> addWidget(add_camera_light_button);

    left_dockwidget_container_layout -> addWidget(light_creation_group);
    left_dockwidget_container -> setLayout(left_dockwidget_container_layout);
    left_dockwidget_container_layout -> addStretch(1);
    this -> left_dockwidget -> setWidget(left_dockwidget_container);
    this -> addDockWidget(Qt::LeftDockWidgetArea, this -> left_dockwidget);
    this -> left_dockwidget -> hide();

    connect(add_scene_light_button,SIGNAL(clicked()),this,SLOT(add_scene_light_slot()));
    connect(add_head_light_button,SIGNAL(clicked()),this,SLOT(add_head_light_slot()));
    connect(add_camera_light_button,SIGNAL(clicked()),this,SLOT(add_camera_light_slot()));


}



void Mainwindow::add_scene_light_slot(){
    this -> add_light(0);
}
void Mainwindow::add_head_light_slot(){
    this -> add_light(1);

}
void Mainwindow::add_camera_light_slot(){
    this -> add_light(2);

}









void Mainwindow::set_action_status(bool enabled, QAction * action) {
    action -> setEnabled(enabled);
}

void Mainwindow::update_GUI_changed_prop() {

    // The status bar is updated depending on whether a shape model remains
    if (this -> wrapped_shape_data.size() == 0 && this -> wrapped_trajectory_data.size() == 0 && this -> wrapped_spacecraft_data.size() == 0) {
        this -> statusBar() -> showMessage("Ready");
    }

    else {

        int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
        std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

        if (this -> wrapped_shape_data.find(name) != this -> wrapped_shape_data.end()){

            auto active_shape_polydata  =  this -> wrapped_shape_data[name] -> get_polydata();
            auto N_vertices = active_shape_polydata -> GetPoints() -> GetNumberOfPoints();
            auto N_facets = active_shape_polydata -> GetNumberOfCells();
            std::string message("Facets : " + std::to_string(N_facets) + " Vertices: " + std::to_string(N_vertices));


            this -> statusBar() -> showMessage(QString::fromStdString(message));
        }
        else if (this -> wrapped_trajectory_data.find(name) != this -> wrapped_trajectory_data.end()){

            auto points =  this -> wrapped_trajectory_data[name] -> get_points();
            std::string message("Trajectory points : " + std::to_string(points -> GetNumberOfPoints()) );

            this -> statusBar() -> showMessage(QString::fromStdString(message));
        }

        else if (this -> wrapped_spacecraft_data.find(name) != this -> wrapped_spacecraft_data.end()){

            std::string message("Spacecraft model: " + name);

            this -> statusBar() -> showMessage(QString::fromStdString(message));
        }

    }

    this -> update_actions_availability();

}



void Mainwindow::open_settings_window() {

    SettingsWindow settings_window(this);
    settings_window.exec();

}

void Mainwindow::createActions() {


    this -> add_small_body_action = new QAction(tr("Load Shape Model"), this);
    this -> add_small_body_action -> setStatusTip(tr("Load obj file holding the facet/vertex description of a shape of interest"));
    connect(this -> add_small_body_action, &QAction::triggered, this, &Mainwindow::add_small_body);

    this -> add_trajectory_action = new QAction(tr("Load Trajectory"), this);
    this -> add_trajectory_action -> setStatusTip(tr("Load a text file storing the x/y/z components a body-fixed trajectory "));
    connect(this -> add_trajectory_action, &QAction::triggered, this, &Mainwindow::add_trajectory);


    this -> open_settings_window_action = new QAction(tr("Settings"), this);
    this -> open_settings_window_action -> setStatusTip(tr("Open settings window where SbgatGUI settings can be set"));
    connect(this -> open_settings_window_action, &QAction::triggered, this, &Mainwindow::open_settings_window);

    this -> show_right_dockwidget_action = new QAction(tr("Show Right Prop Widget"), this);
    this -> show_right_dockwidget_action -> setStatusTip(tr("Shows/hides the right lateral widget holding prop information"));
    connect(this -> show_right_dockwidget_action, &QAction::triggered, this, &Mainwindow::show_right_dockwidget);

    this -> show_left_dockwidget_action = new QAction(tr("Show Left Tool Widget"), this);
    this -> show_left_dockwidget_action -> setStatusTip(tr("Shows/hides lateral widget holding tools menus"));
    connect(this ->show_left_dockwidget_action, &QAction::triggered, this, &Mainwindow::show_left_dockwidget);


    this -> clear_console_action = new QAction(tr("Clear log console"), this);
    this -> clear_console_action -> setStatusTip(tr("Clears the log console"));
    connect(this -> clear_console_action, &QAction::triggered, this, &Mainwindow::clear_console);


    this -> save_console_action = new QAction(tr("Save Log Console"), this);
    this -> save_console_action -> setStatusTip(tr("Saves log console to a file"));
    connect(this -> save_console_action, &QAction::triggered, this, &Mainwindow::save_console);



    this -> compute_geometric_measures_action = new QAction(tr("Compute Geometric Measures"), this);
    this -> compute_geometric_measures_action -> setStatusTip(tr("Compute geometric measures of the selected prop to the console"));
    connect(this -> compute_geometric_measures_action, &QAction::triggered, this, &Mainwindow::compute_geometric_measures);


    this -> compute_pgm_acceleration_action = new QAction(tr("Compute PGM Acceleration"), this);
    this -> compute_pgm_acceleration_action -> setStatusTip(tr("Compute PGM acceleration at a point whose coordinates are expressed in the shape's body frame"));
    connect(this -> compute_pgm_acceleration_action, &QAction::triggered, this, &Mainwindow::compute_pgm_acceleration);


    this -> compute_global_pgm_acceleration_action = new QAction(tr("Compute Global PGM Accelerations"), this);
    this -> compute_global_pgm_acceleration_action -> setStatusTip(tr("Compute PGM accelerations over the entire shape model"));
    connect(this -> compute_global_pgm_acceleration_action, &QAction::triggered, this, &Mainwindow::compute_global_pgm_acceleration);


    this -> compute_global_pgm_potential_action = new QAction(tr("Compute Global PGM Potentials"), this);
    this -> compute_global_pgm_potential_action -> setStatusTip(tr("Compute PGM potentials over the entire shape model"));
    connect(this -> compute_global_pgm_potential_action, &QAction::triggered, this, &Mainwindow::compute_global_pgm_potential);

    this -> compute_grav_slopes_action = new QAction(tr("Compute Gravity Slopes"), this);
    this -> compute_grav_slopes_action -> setStatusTip(tr("Compute PGM accelerations over the entire shape model"));
    connect(this -> compute_grav_slopes_action, &QAction::triggered, this, &Mainwindow::compute_gravity_slopes);

    this -> show_grav_slopes_action = new QAction(tr("Show Gravity Slopes"), this);
    this -> show_grav_slopes_action -> setStatusTip(tr("Display gravity slopes along with colorbar"));
    connect(this -> show_grav_slopes_action, &QAction::triggered, this, &Mainwindow::show_grav_slopes);


    this -> show_global_pgm_pot_action = new QAction(tr("Show Gravity Potentials"), this);
    this -> show_global_pgm_pot_action -> setStatusTip(tr("Display gravity potentials along with colorbar"));
    connect(this -> show_global_pgm_pot_action, &QAction::triggered, this, &Mainwindow::show_global_pgm_pot);

    this -> add_spacecraft_action= new QAction(tr("Load Spacecraft"), this);
    this -> add_spacecraft_action -> setStatusTip(tr("Load spacecraft shape model"));
    connect(this -> add_spacecraft_action, &QAction::triggered, this, &Mainwindow::add_spacecraft);


    this -> move_along_traj_action= new QAction(tr("Move Spacecraft"), this);
    this -> move_along_traj_action -> setStatusTip(tr("Move spacecraft along trajectory"));
    connect(this -> move_along_traj_action, &QAction::triggered, this, &Mainwindow::open_move_along_traj_window);

    this -> open_rendering_properties_window_action =new QAction(tr("Rendering Properties"), this);
    this -> open_rendering_properties_window_action -> setStatusTip(tr("Open window enabling one to change the rendering properties"));
    connect(this -> open_rendering_properties_window_action, &QAction::triggered, this, &Mainwindow::open_rendering_properties_window);

    this -> open_alignment_window_action = new QAction(tr("Actor Alignment and Positioning"),this);
    this -> open_alignment_window_action -> setStatusTip(tr("Open window enabling one rotate or translated any of the displayed props to the center of mass/principal axes of underlying models   "));
    connect(this -> open_alignment_window_action, &QAction::triggered, this, &Mainwindow::open_alignment_window);
    
}

void Mainwindow::update_actions_availability() {

    if (this -> wrapped_shape_data.size() == 0 && this -> wrapped_trajectory_data.size() == 0){
        this -> compute_geometric_measures_action -> setEnabled(false);
        
    }
    else{
        this -> compute_geometric_measures_action -> setEnabled(true);
    }


    if (this -> wrapped_trajectory_data.size() == 0 || this -> wrapped_spacecraft_data.size()== 0){
        this -> move_along_traj_action -> setEnabled(false);
    }

    else if (this -> wrapped_trajectory_data.size()  != 0 && this -> wrapped_spacecraft_data.size() != 0){
        this -> move_along_traj_action -> setEnabled(true);
    }


    if (this -> wrapped_shape_data.size() == 0) {

        this -> compute_pgm_acceleration_action -> setEnabled(false);
        this -> compute_global_pgm_potential_action -> setEnabled(false);
        this -> compute_global_pgm_acceleration_action -> setEnabled(false);
        this -> compute_grav_slopes_action -> setEnabled(false);
    }
    else {

        int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
        std::string name = this -> prop_table -> item(selected_row_index, 0)-> text() .toStdString();

        if (this -> wrapped_shape_data.find(name) != this -> wrapped_shape_data.end()){

            this -> compute_pgm_acceleration_action -> setEnabled(true);
            this -> compute_global_pgm_acceleration_action -> setEnabled(true);
            this -> compute_global_pgm_potential_action -> setEnabled(true);

            if (this -> wrapped_shape_data[name] -> get_global_pgm_acc() == true) {
                this -> compute_grav_slopes_action -> setEnabled(true);
            }

            else {
                this -> compute_grav_slopes_action -> setEnabled(false);
            }
        }
        else{

         this -> compute_pgm_acceleration_action -> setEnabled(false);
         this -> compute_global_pgm_acceleration_action -> setEnabled(false);
         this -> compute_global_pgm_potential_action -> setEnabled(false);

     }
 }






}

void Mainwindow::remove_results_visual_props(std::string name, bool remove_all) {

    // Loop over all the shape mappers to turn off the scalar visibility
    for (auto it = this -> wrapped_shape_data.begin(); it != this -> wrapped_shape_data.end(); ++it) {

        if ((it -> first == name && it -> second -> get_mapper() -> GetScalarVisibility()) || remove_all) {
            it -> second -> get_mapper() -> ScalarVisibilityOff();
            if (this -> renderer -> GetActors2D() -> GetNumberOfItems () > 0) {
                this -> renderer -> RemoveActor2D(this -> renderer -> GetActors2D() -> GetLastActor2D());
            }

            break;
        }

    }


    this -> qvtkWidget -> GetRenderWindow() -> Render();


}

void Mainwindow::show_grav_slopes() {

    QStringList valid_shapes;

    for (auto it = this -> wrapped_shape_data.begin(); it != this -> wrapped_shape_data.end(); ++it) {
        if (it -> second -> get_grav_slopes() == true) {
            valid_shapes << QString::fromStdString(it -> first);
        }
    }

    valid_shapes << QString::fromStdString("");

    bool ok_item;
    QString selected_shape_model = QInputDialog::getItem(this, "Gravitational slopes", "Toggle visibility of gravity slopes for shape model:", valid_shapes, 0, false, &ok_item);


    // If the selected shape model is a valid one
    if (ok_item && this -> wrapped_shape_data.find(selected_shape_model.toStdString()) !=  this -> wrapped_shape_data.end()) {

        // All props are just removed before displaying the one that was just selected
        this -> remove_results_visual_props("", true);

        auto active_mapper  =  this -> wrapped_shape_data[selected_shape_model.toStdString()] -> get_mapper();
        auto active_polydata  =  this -> wrapped_shape_data[selected_shape_model.toStdString()] -> get_polydata();

        if (!active_mapper -> GetScalarVisibility()) {
            active_mapper -> ScalarVisibilityOn();
            active_mapper -> SetScalarModeToUseCellData();

            double range[2] ;
            active_polydata -> GetCellData() -> SetActiveScalars("SlopeData");

            active_polydata -> GetCellData() -> GetScalars() -> GetRange(range);
            active_mapper -> SetColorModeToMapScalars();
            active_mapper -> SetScalarRange(range[0], range[1]);
            vtkLookupTable * lut = vtkLookupTable::SafeDownCast ( active_mapper -> GetLookupTable());
            lut -> SetHueRange(.667, 0);

            vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
            scalarBar -> SetUnconstrainedFontSize (true);
            scalarBar -> GetTitleTextProperty() -> SetFontSize(30);
            scalarBar -> GetLabelTextProperty() -> SetFontSize(30);
            scalarBar -> SetLookupTable(active_mapper -> GetLookupTable());
            scalarBar -> SetTitle("Gravity slopes (deg)");
            scalarBar -> SetNumberOfLabels(4);

            this -> renderer -> AddActor2D(scalarBar);
        }

    }




    this -> qvtkWidget -> GetRenderWindow() -> Render();

}


void Mainwindow::show_global_pgm_pot() {

    QStringList valid_shapes;

    for (auto it = this -> wrapped_shape_data.begin(); it != this -> wrapped_shape_data.end(); ++it) {
        if (it -> second -> get_global_pgm_pot() == true) {
            valid_shapes << QString::fromStdString(it -> first);
        }
    }

    valid_shapes << QString::fromStdString("");


    bool ok_item;
    QString selected_shape_model = QInputDialog::getItem(this, "Gravitational potentials", "Toggle visibility of gravity potentials for shape model:", valid_shapes, 0, false, &ok_item);


    // If the selected shape model is a valid one
    if (ok_item && this -> wrapped_shape_data.find(selected_shape_model.toStdString()) !=  this -> wrapped_shape_data.end()) {

        // All props are just removed before displaying the one that was just selected
        this -> remove_results_visual_props("", true);

        auto active_mapper  =  this -> wrapped_shape_data[selected_shape_model.toStdString()] -> get_mapper();
        auto active_polydata  =  this -> wrapped_shape_data[selected_shape_model.toStdString()] -> get_polydata();

        if (!active_mapper -> GetScalarVisibility()) {
            active_mapper -> ScalarVisibilityOn();
            active_mapper -> SetScalarModeToUseCellData();

            double range[2] ;
            active_polydata -> GetCellData() -> SetActiveScalars("PotentialData");
            active_polydata -> GetCellData() -> GetScalars() -> GetRange(range);
            active_mapper -> SetColorModeToMapScalars();
            active_mapper -> SetScalarRange(range[0], range[1]);
            vtkLookupTable * lut = vtkLookupTable::SafeDownCast ( active_mapper -> GetLookupTable());
            lut -> SetHueRange(.667, 0);



            vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
            scalarBar -> SetUnconstrainedFontSize (true);
            scalarBar -> GetTitleTextProperty() -> SetFontSize(30);
            scalarBar -> GetLabelTextProperty() -> SetFontSize(30);
            scalarBar -> SetLookupTable(active_mapper -> GetLookupTable());
            scalarBar -> SetTitle("Gravity potentials (J)");
            scalarBar -> SetNumberOfLabels(4);

            this -> renderer -> AddActor2D(scalarBar);
        }

    }

    this -> qvtkWidget -> GetRenderWindow() -> Render();
}

void Mainwindow::clear_console() {
    this -> log_console -> clear();
}


void Mainwindow::save_console() {

    QString fileName = QFileDialog::getSaveFileName(this, tr("Save to file"), "",
     tr("Text file (*.txt)"));
    if (fileName != "") {
        QFile file(fileName);

        if (file.open(QIODevice::ReadWrite)) {
            QTextStream stream(&file);
            stream << this -> log_console -> toPlainText();
            file.flush();
            file.close();
        }
        else {
            QMessageBox::critical(this, tr("Error"), tr("Cannot save the file"));
            return;
        }
    }

}


void Mainwindow::show_right_dockwidget() {

    if (this -> right_dockwidget -> isVisible()) {
        this -> right_dockwidget -> hide();
        this -> show_right_dockwidget_action -> setText(QString::fromStdString("Show Prop Widget"));
    }
    else {
        this -> right_dockwidget -> show();
        this -> show_right_dockwidget_action -> setText(QString::fromStdString("Hide Prop Widget"));
    }
}


void Mainwindow::show_left_dockwidget() {

    if (this -> left_dockwidget -> isVisible()) {
        this -> left_dockwidget -> hide();
        this -> show_left_dockwidget_action -> setText(QString::fromStdString("Show Tools Widget"));
    }
    else {
        this -> left_dockwidget -> show();
        this -> show_left_dockwidget_action -> setText(QString::fromStdString("Hide Tools Widget"));
    }
}

void Mainwindow::add_small_body() {

    QString fileName = QFileDialog::getOpenFileName(this,tr("Load Small Body Shape Model"), "~/", tr("Wavefront file (*.obj)"));

    if (fileName.isEmpty() == false) {

        bool ok;
        double scaling_factor = QInputDialog::getDouble(this,"Scaling factor", "Enter scaling factor :", 1, 1e-6, 1e6,5, &ok);
        if (ok) {



         std::stringstream ss;
         ss.str(std::string());

         std::string opening_line = "### Loading shape model ###";
         this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
         this -> log_console -> appendPlainText(QString::fromStdString("- Loading shape model from ") + fileName);

         std::chrono::time_point<std::chrono::system_clock> start, end;

         start = std::chrono::system_clock::now();

            // The name of the shape model is extracted from the path
         int dot_index = fileName.lastIndexOf(".");
         int slash_index = fileName.lastIndexOf("/");
         std::string name = (fileName.toStdString()).substr(slash_index + 1 , dot_index - slash_index - 1);

            // A new ModelDataWrapper is created and stored under the name of the shape model
         std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();

            // The camera is moved to be adjusted to the new shape
         this -> renderer -> GetActiveCamera() -> SetPosition(0, 0, 1.5 * scaling_factor);

            // Reading
         vtkNew<vtkOBJReader> reader;
         reader -> SetFileName(fileName.toStdString().c_str());
         reader -> Update(); 

            // Scaling
         vtkSmartPointer<vtkTransform> transform =
         vtkSmartPointer<vtkTransform>::New();
         transform->Scale(scaling_factor,scaling_factor,scaling_factor);


         vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
         vtkSmartPointer<vtkTransformPolyDataFilter>::New();
         transformFilter->SetInputConnection(reader -> GetOutputPort());
         transformFilter->SetTransform(transform);
         transformFilter -> Update();

            // Create a PolyData
         vtkSmartPointer<vtkPolyData> polygonPolyData = transformFilter -> GetOutput();

            // Create a mapper and actor
         vtkSmartPointer<vtkPolyDataMapper> mapper =
         vtkSmartPointer<vtkPolyDataMapper>::New();

         mapper -> SetInputConnection(reader -> GetOutputPort());
         mapper -> ScalarVisibilityOff();

         vtkSmartPointer<vtkActor> actor =
         vtkSmartPointer<vtkActor>::New();
         actor -> SetMapper(mapper);

            // Visualize
         this -> renderer -> AddActor(actor);

            // Render
         this -> qvtkWidget -> GetRenderWindow() -> Render();

             // Store
         model_data -> set_polydata(polygonPolyData);
         model_data -> set_actor(actor);
         model_data -> set_mapper(mapper);

            // The ModelDataWrapper pointer is stored. 
            // If the name is not already taken, nothing special
         unsigned int count = this -> wrapped_shape_data.count(name);
         if(count == 0){
            this -> wrapped_shape_data[name] = model_data;
        }
            // otherwise, a suffix is added
        else{
            std::string suffix = "(" + std::to_string(count) + ")";
            name = name + suffix;
            this -> wrapped_shape_data[name] = model_data;
        }

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        // The shape table is updated to show the newly loaded shape model
        this -> add_prop_to_table_widget(name);

            // The lateral dockwidget is shown if it was not visible already
        if (this -> right_dockwidget -> isVisible() == false) {
            this -> show_right_dockwidget();

        }

            // The lateral dockwidget is shown if it was not visible already
        if (this -> left_dockwidget -> isVisible() == false) {
            this -> show_left_dockwidget();

        }

            // The GUI actions are updated
        this -> update_actions_availability();

            // The log console displays the name and content of the loaded shape model
        this -> log_console -> appendPlainText(QString::fromStdString("- Loading completed in ")
         + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

        std::string closing_line(opening_line.length() - 1, '#');
        closing_line.append("\n");
        this -> log_console -> appendPlainText(QString::fromStdString(closing_line));

            // A signal is emitted to warn potential dependents
        emit prop_added_signal();

    }
}
}


void Mainwindow::add_spacecraft() {

    QString fileName = QFileDialog::getOpenFileName(this,tr("Load shape model"), "~/", tr("Wavefront file (*.obj)"));

    if (fileName.isEmpty() == false) {

        bool ok;
        double scaling_factor = QInputDialog::getDouble(this,"Scaling factor", "Enter scaling factor :", 1, 1e-6, 1e6,5, &ok);
        if (ok) {


         std::stringstream ss;
         ss.str(std::string());

         std::string opening_line = "### Loading spacecraft model ###";
         this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
         this -> log_console -> appendPlainText(QString::fromStdString("- Loading spacecraft model from ") + fileName);


         std::chrono::time_point<std::chrono::system_clock> start, end;

         start = std::chrono::system_clock::now();


            // The name of the shape model is extracted from the path
         int dot_index = fileName.lastIndexOf(".");
         int slash_index = fileName.lastIndexOf("/");
         std::string name = (fileName.toStdString()).substr(slash_index + 1 , dot_index - slash_index - 1);

         vtkSmartPointer<vtkOBJReader> reader =
         vtkSmartPointer<vtkOBJReader>::New();
         reader -> SetFileName(fileName.toStdString().c_str());
         reader -> Update();

         vtkSmartPointer<vtkPolyData> spacecraft_polydata;



            // The spacececraft actor is aligned with 
           // the center of mass of the loaded spacecraft shape model

         vtkSmartPointer<SBGATMassProperties> center_of_mass_filter =
         vtkSmartPointer<SBGATMassProperties>::New();

         center_of_mass_filter -> SetInputConnection(reader -> GetOutputPort());
         center_of_mass_filter -> Update();

         double center[3];
         center_of_mass_filter -> GetCenterOfMass(center);

         vtkSmartPointer<vtkTransform> translation = vtkSmartPointer<vtkTransform>::New();
         translation -> Translate( - center[0],  - center[1],  - center[2]);
         vtkSmartPointer<vtkTransformPolyDataFilter> translation_filter =
         vtkSmartPointer<vtkTransformPolyDataFilter>::New();
         translation_filter -> SetInputConnection(reader -> GetOutputPort());
         translation_filter -> SetTransform(translation);
         translation_filter -> Update();

         vtkSmartPointer<vtkTransform> scaling = vtkSmartPointer<vtkTransform>::New();
         scaling -> Scale( scaling_factor,  scaling_factor,  scaling_factor);
         vtkSmartPointer<vtkTransformPolyDataFilter> scaling_filter =
         vtkSmartPointer<vtkTransformPolyDataFilter>::New();
         scaling_filter -> SetInputConnection(translation_filter -> GetOutputPort());
         scaling_filter -> SetTransform(scaling);
         scaling_filter -> Update();




            // Set the actor and mapper
         vtkSmartPointer<vtkPolyDataMapper> mapper = 
         vtkSmartPointer<vtkPolyDataMapper>::New();
         mapper -> SetInputConnection(scaling_filter->GetOutputPort());
         vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
         actor -> SetMapper(mapper);
         this -> renderer -> AddActor(actor);

            // A new ModelDataWrapper is created and stored under the name of the spacecraft
         std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();
         model_data -> set_polydata(mapper -> GetInput());
         model_data -> set_actor(actor);
         model_data -> set_mapper(mapper);


            // The ModelDataWrapper pointer is stored. 
            // If the name is not already taken, nothing special
         unsigned int count = this -> wrapped_spacecraft_data.count(name);
         if(count == 0){
            this -> wrapped_spacecraft_data[name] = model_data;
        }
            // otherwise, a suffix is added
        else{
            std::string suffix = "(" + std::to_string(count) + ")";
            name = name + suffix;
            this -> wrapped_spacecraft_data[name] = model_data;
        }



        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

            // The shape table is updated to show the newly loaded shape model
        this -> add_prop_to_table_widget(name);


            // The lateral dockwidget is shown if it was not visible already
        if (this -> right_dockwidget -> isVisible() == false) {
            this -> show_right_dockwidget();

        }

            // The GUI actions are updated
        this -> update_actions_availability();

                  // The log console displays the name and content of the loaded shape model
        this -> log_console -> appendPlainText(QString::fromStdString("- Loading completed in ")
         + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

        std::string closing_line(opening_line.length() - 1, '#');
        closing_line.append("\n");
        this -> log_console -> appendPlainText(QString::fromStdString(closing_line));


            // A signal is emitted to warn potential dependents
        emit prop_added_signal();

    }
}
}



void Mainwindow::add_trajectory() {

    QString fileName = QFileDialog::getOpenFileName(this,
     tr("Load trajectory"), "~/", tr("(*.txt)"));

    if (fileName.isEmpty() == false) {

        bool ok;
        double scaling_factor = QInputDialog::getDouble(this,
            "Scaling factor", "Enter scaling factor :", 1, 1e-6, 1e6,
            5, &ok);
        if (ok) {


         std::stringstream ss;
         ss.str(std::string());

         std::string opening_line = "### Loading spacecraft model ###";
         this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
         this -> log_console -> appendPlainText(QString::fromStdString("- Loading trajectory from ") + fileName);

         std::chrono::time_point<std::chrono::system_clock> start, end;

         start = std::chrono::system_clock::now();


            // The name of the trajectory is extracted from the path
         int dot_index = fileName.lastIndexOf(".");
         int slash_index = fileName.lastIndexOf("/");
         std::string name = (fileName.toStdString()).substr(slash_index + 1 , dot_index - slash_index - 1);

         arma::mat traj;
         traj.load(fileName.toStdString());

         traj *= scaling_factor;

         if(traj.n_cols < traj.n_rows){
            arma::inplace_trans(traj);
        }


            // Create a vtkPoints object and store the points in it
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

        for (unsigned int i = 0; i < traj.n_cols; ++i){

            arma::vec pos = traj.col(i).rows(1,3);
            points->InsertNextPoint(pos.colptr(0));

        }

        vtkSmartPointer<vtkParametricSpline> spline = 
        vtkSmartPointer<vtkParametricSpline>::New();
        spline -> SetPoints(points);

        vtkSmartPointer<vtkParametricFunctionSource> functionSource = 
        vtkSmartPointer<vtkParametricFunctionSource>::New();
        functionSource -> SetParametricFunction(spline);
        functionSource -> SetUResolution (traj.n_cols);
        functionSource -> Update();

            // Set the actor and mapper
        vtkSmartPointer<vtkPolyDataMapper> mapper = 
        vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper -> SetInputConnection(functionSource -> GetOutputPort());
        vtkSmartPointer<vtkActor> actor = 
        vtkSmartPointer<vtkActor>::New();
        actor -> SetMapper(mapper);
        this -> renderer -> AddActor(actor);

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

            // A new ModelDataWrapper is created and stored under the name of the trajectory
        std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();
        model_data -> set_polydata(mapper -> GetInput());
        model_data -> set_points(points);
        model_data -> set_actor(actor);
        model_data -> set_mapper(mapper);


            // The ModelDataWrapper pointer is stored. 
            // If the name is not already taken, nothing special
        unsigned int count = this -> wrapped_trajectory_data.count(name);
        if(count == 0){
            this -> wrapped_trajectory_data[name] = model_data;
        }
            // otherwise, a suffix is added
        else{
            std::string suffix = "(" + std::to_string(count) + ")";
            name = name + suffix;
            this -> wrapped_trajectory_data[name] = model_data;
        }


            // The shape table is updated to show the newly loaded shape model
        this -> add_prop_to_table_widget(name);

            // The lateral dockwidget is shown if it was not visible already
        if (this -> right_dockwidget -> isVisible() == false) {
            this -> show_right_dockwidget();
        }

            // The GUI actions are updated
        this -> update_actions_availability();

            // The log console displays the name and content of the loaded shape model
        this -> log_console -> appendPlainText(QString::fromStdString("- Loading completed in ")
         + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

            // A signal is emitted to warn potential dependents
        emit prop_added_signal();

    }
}
}






void Mainwindow::add_prop_to_table_widget(std::string name) {


    QTableWidgetItem * nameItem = new QTableWidgetItem(QString::fromStdString(name));
    nameItem -> setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);


    this -> prop_table -> insertRow(this -> prop_table -> rowCount());
    this -> prop_table -> setItem(this -> prop_table -> rowCount() - 1, 0, nameItem);


    QTableWidgetItem *checkBoxItem = new QTableWidgetItem();
    checkBoxItem -> setCheckState(Qt::Checked);
    
    this -> prop_table -> setItem(this -> prop_table -> rowCount() - 1, 1, checkBoxItem);

    QWidget * button_container = new QWidget(this -> prop_table -> cellWidget(this -> prop_table -> rowCount() - 1, 0));
    QHBoxLayout* layout = new QHBoxLayout(button_container);

    QPushButton * erase_shape_button = new QPushButton(button_container);
    erase_shape_button -> setText("X");
    
    erase_shape_button -> setProperty("name", QVariant(QString::fromStdString(name)));
    layout -> addWidget(erase_shape_button);
    layout -> setAlignment(Qt::AlignCenter);
    layout -> setContentsMargins(0, 0, 0, 0);
    button_container -> setLayout(layout);

    this -> prop_table -> setCellWidget(this -> prop_table -> rowCount() - 1, 3, button_container);

    // The push button is connected to the proper slot
    connect(erase_shape_button, SIGNAL(clicked(bool)), this, SLOT(remove_prop()));

    // The check box is connected to the proper slot
    connect(this -> prop_table, SIGNAL(cellChanged(int, int)), this, SLOT(toggle_prop_visibility(int, int)));

    // This prop is selected
    this -> prop_table -> selectRow(this -> prop_table -> rowCount() - 1);


}


void Mainwindow::toggle_prop_visibility(int row, int col) {


    if (col == 1) {
        std::string name = this -> prop_table -> item(row, 0) -> text() .toStdString();

        auto item = this -> prop_table -> item(row, col);


        if ( this -> wrapped_shape_data.find(name)!= this -> wrapped_shape_data.end()){
            if (item -> checkState() == Qt::Checked) {
                this -> wrapped_shape_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_actor() -> VisibilityOn();
            }

            else {
                this -> wrapped_shape_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_actor() -> VisibilityOff();
                this -> remove_results_visual_props(this -> prop_table -> item(row, 0) -> text() . toStdString(), false);
            }
        }
        else if ( this -> wrapped_trajectory_data.find(name) != this -> wrapped_trajectory_data.end()){
         if (item -> checkState() == Qt::Checked) {
            this -> wrapped_trajectory_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_actor() -> VisibilityOn();
        }

        else {
            this -> wrapped_trajectory_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_actor() -> VisibilityOff();
            this -> remove_results_visual_props(this -> prop_table -> item(row, 0) -> text() . toStdString(), false);
        }
    }
    else if ( this -> wrapped_light_data.find(name) != this -> wrapped_light_data.end()){
     if (item -> checkState() == Qt::Checked) {
        this -> wrapped_light_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_light_actor() -> SetVisibility(1);
    }

    else {
        this -> wrapped_light_data[this -> prop_table -> item(row, 0) -> text() . toStdString()] -> get_light_actor() -> SetVisibility(0);
        this -> remove_results_visual_props(this -> prop_table -> item(row, 0) -> text() . toStdString(), false);
    }
}



}

    // The Render window is updated
this -> qvtkWidget -> GetRenderWindow() -> Render();


}


void Mainwindow::remove_prop() {

    QPushButton * button = qobject_cast<QPushButton*>(sender());
    std::string name = button -> property("name") . toString().toStdString();


    if (this -> wrapped_shape_data.find(name) != this -> wrapped_shape_data.end() ){
        // The actor of this shape is removed
        this -> renderer -> RemoveActor(this -> wrapped_shape_data[name] -> get_actor());

        // The data wrapper is removed
        this -> wrapped_shape_data.erase(name);
    }
    else if (this -> wrapped_trajectory_data.find(name) != this -> wrapped_trajectory_data.end() ){
        // The actor of this trajectory is removed
        this -> renderer -> RemoveActor(this -> wrapped_trajectory_data[name] -> get_actor());

        // The data wrapper is removed
        this -> wrapped_trajectory_data.erase(name);
    }
    else if (this -> wrapped_spacecraft_data.find(name) != this -> wrapped_spacecraft_data.end()){
        this -> renderer -> RemoveActor(this -> wrapped_spacecraft_data[name] -> get_actor());

        // The data wrapper is removed
        this -> wrapped_spacecraft_data.erase(name);
    }
    else if (this -> wrapped_light_data.find(name) != this -> wrapped_light_data.end()){

        if (this -> wrapped_light_data.size() < 2){
            // If there's only one light left, it cannot be removed as VTK will not 
            // allow this
            return;
        }

        // The light actor is removed
        this -> renderer -> RemoveActor(this -> wrapped_light_data[name] -> get_light_actor());

        // The light is removed from the renderer
        this -> get_renderer() -> RemoveLight(this -> wrapped_light_data[name] -> get_light());

        // The data wrapper is removed
        this -> wrapped_light_data.erase(name);


    }

    // The corresponding row in the table widget is removed
    // This will trigger the corresponding signal/slot mechanism updating the GUI
    for (int i = 0; i < this -> prop_table -> rowCount(); ++i) {
        if (this -> prop_table -> item(i, 0) -> text() == QString::fromStdString(name)) {
            this -> prop_table -> removeRow(i);
            break;
        }
    }

    this -> update_actions_availability();

    // The Render window is updated
    this -> qvtkWidget -> GetRenderWindow() -> Render();

    emit prop_removed_signal();

}



void Mainwindow::compute_geometric_measures(){

 int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
 std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

 if (this -> wrapped_shape_data.find(name) != this -> wrapped_shape_data.end()){

    std::stringstream ss;

    ss.str(std::string());
    ss.precision(10);
    
    std::string opening_line = "### Computing shape model geometric measures ###";
    this -> log_console -> appendPlainText(QString::fromStdString(opening_line));

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    vtkSmartPointer<SBGATMassProperties> mass_properties_filter = vtkSmartPointer<SBGATMassProperties>::New();
    mass_properties_filter -> SetInputData(this -> wrapped_shape_data[name] -> get_polydata());
    mass_properties_filter -> Update();
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;


    this -> log_console -> appendPlainText(QString::fromStdString("\n- Surface of " + name + " (m^2) :"));
    this -> log_console -> appendPlainText(" " + QString::number(mass_properties_filter -> GetSurfaceArea ()));

    this -> log_console -> appendPlainText(QString::fromStdString("\n- Volume of " + name + " (m^3) :"));
    this -> log_console -> appendPlainText(" " + QString::number(mass_properties_filter -> GetVolume()));

    this -> log_console -> appendPlainText(QString::fromStdString("\n- Bounding box of " + name + " (m) :"));

    double * bbox =  mass_properties_filter -> GetBoundingBox();

    this -> log_console -> appendPlainText(QString::fromStdString("-- Min: " + std::to_string(bbox[0]) + " "+ std::to_string(bbox[2]) + " "+ std::to_string(bbox[4])));
    this -> log_console -> appendPlainText(QString::fromStdString("-- Max: " + std::to_string(bbox[1]) + " "+ std::to_string(bbox[3]) + " "+ std::to_string(bbox[5])));

    ss.str(std::string());
    ss.precision(10);

    this -> log_console -> appendPlainText(QString::fromStdString("\n- Center of mass of " + name + " (m) :"));
    mass_properties_filter -> GetCenterOfMass().t().raw_print(ss);
    this -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

    ss.str(std::string());
    ss.precision(10);

    this -> log_console -> appendPlainText(QString::fromStdString("- Dimensionless inertia tensor of " + name ));
    mass_properties_filter -> GetInertiaTensor().raw_print(ss);
    this -> log_console -> appendPlainText(QString::fromStdString(ss.str()));
    
    ss.str(std::string());
    ss.precision(10);

    this -> log_console -> appendPlainText(QString::fromStdString("- Principal axes of " + name ));
    mass_properties_filter -> GetPrincipalAxes().raw_print(ss);
    this -> log_console -> appendPlainText(QString::fromStdString(ss.str()));

    ss.str(std::string());
    ss.precision(10);

    this -> log_console -> appendPlainText(QString::fromStdString("- Dimensionless inertia moments of " + name ));
    mass_properties_filter -> GetInertiaMoments().t().raw_print(ss);
    this -> log_console -> appendPlainText(QString::fromStdString(ss.str()));


    this -> log_console -> appendPlainText(QString::fromStdString("- Done computing in ")
     + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

    std::string closing_line(opening_line.length() - 1, '#');
    closing_line.append("\n");
    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));

}
else if (this -> wrapped_trajectory_data.find(name) != this -> wrapped_trajectory_data.end()){
    std::string opening_line = "### Computing trajectory geometry measures ###\n";
    this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
    std::string closing_line(opening_line.length() - 1, '#');
    closing_line.append("\n");
    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
}

}


void Mainwindow::add_light(int light_type){


    std::vector<std::string> light_names = {"Scene light","Head light","Camera light"};
    vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();


    light -> SetFocalPoint(this  -> get_renderer() -> GetActiveCamera() -> GetFocalPoint());
    light -> SetPosition(this  -> get_renderer() -> GetActiveCamera() -> GetPosition());
    light -> SetColor(1,1,1);
    light -> SetIntensity(1);
    light -> SetLightType(light_type); // 0: scene , 1: Head light , 2: camera light

    // The light and its props are stored in a ModelDataWrapper
    std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();
    vtkSmartPointer<vtkLightActor> lightActor = vtkSmartPointer<vtkLightActor>::New();
    lightActor -> SetLight(light);    
    model_data -> set_light(light);
    model_data -> set_light_actor(lightActor);

    light->SetPositional(true); // required for vtkLightActor below
    light->SetConeAngle(20);


    // At this stage, name does not account for repeated lights of the same type
    std::string name = light_names[light_type];

    // A potential suffix to add to the light name is found by 
    // looking at lights that already exist
    int light_count = 0;
    for (auto light_it = this -> wrapped_light_data.begin(); light_it != this -> wrapped_light_data.end(); ++light_it){

        if (light_it -> second -> get_light() -> GetLightType() == light_type){
            ++light_count;
        }
    }

    if(light_count == 0){
        this -> wrapped_light_data[name] = model_data;
    }
            // otherwise, a suffix is added
    else{
        std::string suffix = "(" + std::to_string(light_count) + ")";
        name = name + suffix;
        this -> wrapped_light_data[name] = model_data;
    }

    this -> get_renderer() -> AddViewProp(lightActor);
    this -> get_renderer() -> AddLight(light);

    // The prop table is updated to show the newly loaded prop
    this -> add_prop_to_table_widget(name);

    this -> qvtkWidget -> GetRenderWindow() -> Render();




}


void Mainwindow::compute_global_pgm_acceleration() {

    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    // auto dyn_analyses = std::make_shared<SBGAT_CORE::DynamicAnalyses>(this -> wrapped_shape_data[name] -> get_shape_model().get());

    // bool ok_density = false;


    // double density = QInputDialog::getDouble(this,
    //  "Global Polyhedron Gravity Model Acceleration", "Density (kg/m^3) :", 2000, 0, 1e9,
    //  5, &ok_density);




    // if (ok_density) {

    //     double mu = density * arma::datum::G * this -> wrapped_shape_data[name] -> get_shape_model() -> get_volume();


    //     // Log in
    //     this -> log_console -> appendPlainText(QString::fromStdString("- Computing global PGM facet accelerations of " + name + " ..."));

    //     // The shape selection table is frozen
    //     this -> prop_table -> setDisabled(true);
    //     this -> menuBar() -> setDisabled(true);


    //     QThread * thread = new QThread;
    //     Worker * worker = new Worker(dyn_analyses, mu, this -> wrapped_shape_data[name],
    //      name);
    //     worker ->  moveToThread(thread);
    //     connect(thread, SIGNAL(started()), worker, SLOT(process_pgm_acc()));
    //     connect(worker, SIGNAL(finished()), thread, SLOT(quit()));
    //     connect(worker, SIGNAL(logging_out(QString)), this -> log_console, SLOT(insertPlainText(QString)));
    //     connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    //     connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    //     connect(worker, SIGNAL(finished()), this, SLOT(update_actions_availability()));
    //     connect(worker, SIGNAL(free_shape_table(bool)), this -> prop_table, SLOT(setEnabled(bool)));
    //     connect(worker, SIGNAL(free_menu_bar(bool)), this -> menuBar() , SLOT(setEnabled(bool)));

    //     thread -> start();


    // }

}



void Mainwindow::compute_global_pgm_potential() {

    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    // auto dyn_analyses = std::make_shared<SBGAT_CORE::DynamicAnalyses>(this -> wrapped_shape_data[name] -> get_shape_model().get());

    // bool ok_density = false;


    // double density = QInputDialog::getDouble(this,
    //  "Global Polyhedron Gravity Model Acceleration", "Density (kg/m^3) :", 2000, 0, 1e9,
    //  5, &ok_density);

    // if (ok_density) {


    //     double mu = density * arma::datum::G * this -> wrapped_shape_data[name] -> get_shape_model() -> get_volume();



    //     // Log int
    //     this -> log_console -> appendPlainText(QString::fromStdString("- Computing global PGM facet potentials of " + name + " ..."));

    //     // The shape selection table is frozen
    //     this -> prop_table -> setDisabled(true);
    //     this -> menuBar()  -> setDisabled(true);


    //     QThread * thread = new QThread;
    //     Worker * worker = new Worker(dyn_analyses,
    //      mu,
    //      this -> wrapped_shape_data[name],
    //      name);
    //     worker ->  moveToThread(thread);
    //     connect(thread, SIGNAL(started()), worker, SLOT(process_pgm_pot()));
    //     connect(worker, SIGNAL(finished()), thread, SLOT(quit()));
    //     connect(worker, SIGNAL(logging_out(QString)), this -> log_console, SLOT(insertPlainText(QString)));
    //     connect(worker, SIGNAL(finished()), worker, SLOT(deleteLater()));
    //     connect(thread, SIGNAL(finished()), thread, SLOT(deleteLater()));
    //     connect(worker, SIGNAL(finished()), this, SLOT(update_actions_availability()));
    //     connect(worker, SIGNAL(finished()), this, SLOT(update_vtk_potentials()));
    //     connect(worker, SIGNAL(free_shape_table(bool)), this -> prop_table, SLOT(setEnabled(bool)));
    //     connect(worker, SIGNAL(free_menu_bar(bool)), this -> menuBar() , SLOT(setEnabled(bool)));

    //     thread -> start();

    // }

}


void Mainwindow::compute_gravity_slopes() {
    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();


    // SBGAT_CORE::DynamicAnalyses dynas(this -> wrapped_shape_data[name] -> get_shape_model().get());

    // bool ok_spin_axis = true;
    // bool correct_format = false;
    // bool ok_spin_rate = false;
    // arma::vec spin_axis = {0, 0, 1};
    // arma::vec angles_arma = {0, 0, 0};

    // while (ok_spin_axis == true && correct_format == false) {

    //     QString coords = QInputDialog::getText(this,
    //        tr("Gravity slopes"),
    //        tr("(Azimuth,Elevation) of spin axis (deg) . (0,0) : along z :"),
    //        QLineEdit::Normal,
    //        QString::fromStdString("Azimuth,Elevation") ,
    //        &ok_spin_axis);

    //     QStringList angles_split = coords.split(",");

    //     if (angles_split.count() != 2) {
    //         correct_format = false;
    //         continue;
    //     }

    //     // Matching doubles
    //     QRegExp re("[-+]?[0-9]*\\.?[0-9]+");

    //     if (re.exactMatch(angles_split.at(0)) && re.exactMatch(angles_split.at(1))) {
    //         angles_arma(0) = angles_split.at(0).toDouble();
    //         angles_arma(1) = angles_split.at(1).toDouble();

    //         correct_format = true;


    //     }
    //     else {
    //         correct_format = false;
    //         continue;
    //     }

    // }

    // if (ok_spin_axis) {
    //     spin_axis = RBK::euler313_to_dcm(angles_arma).t() * spin_axis;

    //     double period = QInputDialog::getDouble(this,
    //         "Gravity slopes", "Rotation period (hours) :", 0, -1e9, 1e9,
    //         5, &ok_spin_rate);

    //     double spin_rate = arma::datum::pi * 2 / (period * 3600);

    //     if (ok_spin_rate ) {

    //         this -> log_console -> appendPlainText(QString::fromStdString("- Computing gravity slopes of " + name +  "..."));

    //         std::chrono::time_point<std::chrono::system_clock> start, end;
    //         start = std::chrono::system_clock::now();
    //         arma::vec slope_stats = dynas.compute_gravity_slopes(spin_axis, spin_rate);

    //         this -> update_vtk_slopes();


    //         // A GUI flag is updated to indicate that this shape model has consistent slopes ready to be displayed
    //         this -> wrapped_shape_data[name] -> set_grav_slopes(true);

    //         // The GUI actions are updated
    //         this -> update_actions_availability();

    //         // The previously displayed slopes (if any) are removed
    //         this -> remove_results_visual_props("", true);

    //         // Some statistics regarding the slopes are printed
    //         std::string message("-- Mean slope: " + std::to_string(slope_stats(1)) + " deg");
    //         this -> log_console -> appendPlainText(QString::fromStdString(message));
    //         message = ("-- Minimum slope: " + std::to_string(slope_stats(0)) + " deg");
    //         this -> log_console -> appendPlainText(QString::fromStdString(message));
    //         message = ("-- Maximum slope: " + std::to_string(slope_stats(2)) + " deg");
    //         this -> log_console -> appendPlainText(QString::fromStdString(message));

    //         end = std::chrono::system_clock::now();
    //         std::chrono::duration<double> elapsed_seconds = end - start;

    //         this -> log_console -> appendPlainText(QString::fromStdString("- Done computing in ")
    //            + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

    //     }

    // }

}

void Mainwindow::update_vtk_potentials() {

    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    // vtkSmartPointer<vtkPolyData> active_polydata =  this -> wrapped_shape_data[name] -> get_polydata();
    // vtkSmartPointer<vtkPolyDataMapper> active_mapper = this -> wrapped_shape_data[name] -> get_mapper();
    // std::shared_ptr<SBGAT_CORE::ShapeModel> active_shape_model = this -> wrapped_shape_data[name] -> get_shape_model();

    // vtkSmartPointer<vtkDoubleArray> potential_data =
    // vtkSmartPointer<vtkDoubleArray>::New();
    // potential_data -> SetNumberOfValues(active_shape_model -> get_NFacets());
    // potential_data -> SetName("PotentialData");

    // for (unsigned int facet_index = 0; facet_index < active_shape_model -> get_NFacets(); ++facet_index) {

    //     SBGAT_CORE::Facet * facet =  active_shape_model -> get_facets() -> at(facet_index);

    //     potential_data -> SetValue(facet_index, facet -> get_facet_results() -> get_grav_potential());
    // }

    // // The array is added to the polydata
    // active_polydata -> GetCellData() -> SetActiveScalars("PotentialData");
    // active_polydata -> GetCellData() -> SetScalars(potential_data);
    // active_polydata -> Modified();

}

void Mainwindow::open_move_along_traj_window(){

    MoveAlongTrajectoryWindow * move_along_traj_window = new MoveAlongTrajectoryWindow(this);
    connect(this,SIGNAL(prop_removed_signal()),move_along_traj_window,SLOT(prop_removed_slot()));
    connect(this,SIGNAL(prop_added_signal()),move_along_traj_window,SLOT(prop_added_slot()));
    
    move_along_traj_window -> show();
}   


void Mainwindow::open_rendering_properties_window(){

    RenderingPropertiesWindow * rendering_properties_window = new RenderingPropertiesWindow(this);

    connect(this,SIGNAL(prop_removed_signal()),rendering_properties_window,SLOT(prop_removed_slot()));
    connect(this,SIGNAL(prop_added_signal()),rendering_properties_window,SLOT(prop_added_slot()));
    
    rendering_properties_window -> show();
}   




void Mainwindow::update_vtk_slopes() {

    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    // vtkSmartPointer<vtkPolyData> active_polydata =  this -> wrapped_shape_data[name] -> get_polydata();
    // vtkSmartPointer<vtkPolyDataMapper> active_mapper = this -> wrapped_shape_data[name] -> get_mapper();
    // std::shared_ptr<SBGAT_CORE::ShapeModel> active_shape_model = this -> wrapped_shape_data[name] -> get_shape_model();

    // vtkSmartPointer<vtkDoubleArray> slope_data =
    // vtkSmartPointer<vtkDoubleArray>::New();
    // slope_data -> SetNumberOfValues(active_shape_model -> get_NFacets());
    // slope_data -> SetName("SlopeData");

    // for (unsigned int facet_index = 0; facet_index < active_shape_model -> get_NFacets(); ++facet_index) {

    //     SBGAT_CORE::Facet * facet =  active_shape_model -> get_facets() -> at(facet_index);

    //     slope_data -> SetValue(facet_index, facet -> get_facet_results() -> get_grav_slope());
    // }

    // // The array is added to the polydata
    // active_polydata -> GetCellData() -> SetActiveScalars("SlopeData");
    // active_polydata -> GetCellData() -> SetScalars(slope_data);
    // active_polydata -> Modified();

}


void Mainwindow::compute_pgm_acceleration() {

    // int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    // std::string name = this -> prop_table -> item(selected_row_index, 0)-> text() .toStdString();

    // SBGAT_CORE::DynamicAnalyses dynas(this -> wrapped_shape_data[name] -> get_shape_model().get());

    // bool ok_coords = true;
    // bool correct_format = false;
    // double point[3];
    // bool ok_density = false;
    // double density;

    // while (ok_coords == true && correct_format == false) {

    //     QString coords = QInputDialog::getText(this,
    //        tr("Polyhedron Gravity Model Acceleration"),
    //        tr("Body-fixed frames coordinates (m) :"),
    //        QLineEdit::Normal,
    //        QString::fromStdString("x,y,z") ,
    //        &ok_coords);

    //     QStringList coords_split = coords.split(",");

    //     if (coords_split.count() != 3) {
    //         correct_format = false;
    //         continue;
    //     }

    //     // Matching doubles
    //     QRegExp re("[-+]?[0-9]*\\.?[0-9]+");

    //     if (re.exactMatch(coords_split.at(0)) && re.exactMatch(coords_split.at(1)) && re.exactMatch(coords_split.at(2))) {
    //         point[0] = coords_split.at(0).toDouble();
    //         point[1] = coords_split.at(1).toDouble();
    //         point[2] = coords_split.at(2).toDouble();

    //         correct_format = true;


    //     }
    //     else {
    //         correct_format = false;
    //         continue;
    //     }

    // }



    // if (ok_coords) {
    //     density = QInputDialog::getDouble(this,
    //       "Polyhedron Gravity Model Acceleration", "Density (kg/m^3) :", 2000, 0, 1e9,
    //       5, &ok_density);


    //     if (ok_density) {

    //         double mu = density * arma::datum::G * this -> wrapped_shape_data[name] -> get_shape_model() -> get_volume();



    //         // The PGM acceleration is computed at the provided point
    //         arma::vec coords_arma = {point[0], point[1], point[2]};
    //         std::stringstream ss_coords;
    //         ss_coords.precision(10);
    //         ss_coords << " " << point[0] << "\n" << " " << point[1] << "\n" << " " << point[2] << "\n";

    //         arma::vec acc = dynas.pgm_acceleration(point , mu);
    //         std::stringstream ss_acc;
    //         ss_acc.precision(10);
    //         ss_acc << " " << acc[0] << "\n" << " " << acc[1] << "\n" << " " << acc[2] << "\n";

    //         this -> log_console -> appendPlainText(QString::fromStdString("\n- At body-fixed coordinates (m) : "));
    //         this -> log_console -> appendPlainText(QString::fromStdString(ss_coords.str()));
    //         this -> log_console -> appendPlainText(QString::fromStdString("- PGM acceleration (m/s^2): "));
    //         this -> log_console -> appendPlainText(QString::fromStdString(ss_acc.str()));

    //     }
    // }


}




void Mainwindow::createMenus() {

    this -> SmallBodyMenu = this -> menuBar() -> addMenu(tr("&Small Body"));
    this -> SmallBodyMenu -> addAction(this -> add_small_body_action);
    this -> SmallBodyMenu -> addSeparator();
    this -> SmallBodyMenu -> addAction(this -> open_settings_window_action);


    this -> TrajectoryMenu = this -> menuBar() -> addMenu(tr("&Trajectory"));
    this -> TrajectoryMenu -> addAction(this -> add_trajectory_action);

    this -> SpacecraftMenu = this -> menuBar() -> addMenu(tr("&Spacecraft"));
    this -> SpacecraftMenu -> addAction(this -> add_spacecraft_action);
    this -> SpacecraftMenu -> addAction(this -> move_along_traj_action);

    
    this -> MeasuresMenu = menuBar() -> addMenu(tr("&Measures"));
    this -> MeasuresMenu -> addAction(this -> compute_geometric_measures_action);


    this -> DynamicAnalysesMenu = menuBar() -> addMenu(tr("&Analyses"));
    // this -> DynamicAnalysesMenu -> addAction(this -> compute_pgm_acceleration_action);
    // this -> DynamicAnalysesMenu -> addSeparator();
    // this -> DynamicAnalysesMenu -> addAction(this -> compute_global_pgm_potential_action);
    // this -> DynamicAnalysesMenu -> addAction(this -> compute_global_pgm_acceleration_action);
    // this -> DynamicAnalysesMenu -> addAction(this -> compute_grav_slopes_action);


    this -> ResultsMenu = menuBar() -> addMenu(tr("&Visualization"));
    this -> ResultsMenu -> addAction(this -> open_rendering_properties_window_action);
    this -> ResultsMenu -> addAction(this -> open_alignment_window_action);

    this -> ResultsMenu -> addSeparator();

    this -> ResultsMenu -> addAction(this -> show_grav_slopes_action);
    this -> ResultsMenu -> addAction(this -> show_global_pgm_pot_action);


    this -> ViewMenu = menuBar() -> addMenu(tr("&View"));
    this -> ViewMenu -> addAction(this -> show_left_dockwidget_action);
    this -> ViewMenu -> addAction(this -> show_right_dockwidget_action);
    this -> ViewMenu -> addSeparator();

    this -> ConsoleMenu = menuBar() -> addMenu(tr("&Console"));
    this -> ConsoleMenu -> addAction(this -> clear_console_action);
    this -> ConsoleMenu -> addAction(this -> save_console_action);

}


void Mainwindow::open_alignment_window(){

}


DataMap Mainwindow::get_wrapped_shape_data() const{
    return this -> wrapped_shape_data;
}

DataMap Mainwindow::get_wrapped_trajectory_data() const{
    return this -> wrapped_trajectory_data;
}

DataMap Mainwindow::get_wrapped_attitude_data() const{
    return this -> wrapped_attitude_data;
}

DataMap Mainwindow::get_wrapped_spacecraft_data() const{
    return this -> wrapped_spacecraft_data;
}




std::pair<std::string,vtkSmartPointer<vtkActor> > Mainwindow::get_skybox_pair() const{
    return this ->skybox_pair;
}

void Mainwindow::set_skybox_actor(vtkSmartPointer<vtkActor> skybox_actor){
    this -> skybox_pair.second = skybox_actor;
}


void Mainwindow::set_skybox_directory(std::string skybox_dir){
    this -> skybox_pair.first = skybox_dir;
}



vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}







