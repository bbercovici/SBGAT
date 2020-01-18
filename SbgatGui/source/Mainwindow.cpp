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
#include <vtkPointData.h>

#include <vtkCenterOfMass.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkGenericRenderWindowInteractor.h>
#include <vtkActor2DCollection.h>
#include <vtkCellData.h>
#include <vtkTextProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkProperty.h>
#include <vtkIdFilter.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkPolyDataNormals.h>
#include <SBGATMassProperties.hpp>
#include <SBGATObjWriter.hpp>


#include "SettingsWindow.hpp"
#include "RenderingPropertiesWindow.hpp"
#include "YORPWindow.hpp"
#include "SHARMWindow.hpp"
#include "RadarWindow.hpp"
#include "LCWindow.hpp"
#include "SurfacePGMWindow.hpp"
#include "SelectMapperWindow.hpp"
#include "VertexEditionWindow.hpp"
#include "MassPropertiesWindow.hpp"


using namespace SBGAT_GUI;


Mainwindow::Mainwindow() {

    // The GUI elements are created
    this -> setupUi();

    // The elements interfacing with SbgatCore are created
    this -> frame_graph = std::make_shared<SBGATFrameGraph>();
    this -> frame_graph -> add_frame("inertial_default"); // a default inertial frame of reference is created

}

void Mainwindow::setupUi() {
    // this -> resize(1024, 768);

    // The widget are created
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
    QStringList header_lists = {"Name", "State" , "Visibility", ""};
    this -> prop_table -> setHorizontalHeaderLabels(header_lists);
    this -> prop_table -> horizontalHeader() -> setStretchLastSection(true);
    this -> prop_table -> setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    
    // Selecting an item in the table highlights the entire row
    this -> prop_table -> setSelectionBehavior(QAbstractItemView::SelectRows);
    this -> prop_table -> setSelectionMode(QAbstractItemView::SingleSelection);

    // Prevents edition of labels
    this -> prop_table -> setEditTriggers(QAbstractItemView::NoEditTriggers);

    // Central window
    this -> setCentralWidget(this -> qvtkWidget);
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



    // The lateral dockwidgets are initialized
    // This is delayed until after the renderer is updated so
    // that the default props (default light...) have been instantiated

    this -> show();

    this -> init_right_dockwidget();
    this -> qvtkWidget -> GetRenderWindow() -> Render();




}

void Mainwindow::init_rendering_window(){

 // A VTK renderer is created and linked with the qvtk widget
    this -> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> render_window = vtkSmartPointer< vtkGenericOpenGLRenderWindow>::New();
    this -> qvtkWidget -> SetRenderWindow(render_window);
    this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(this -> renderer);


    this -> renderer -> SetGradientBackground (true);
    this -> renderer -> SetBackground (0.5, 0.5, 1);



    vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
    

    renderWindowInteractor -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow());
    


    vtkSmartPointer<vtkAxesActor> axes = vtkSmartPointer<vtkAxesActor>::New();

    this -> orientation_widget =
    vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    orientation_widget -> SetOrientationMarker( axes );
    orientation_widget -> SetInteractor(render_window -> GetInteractor() );
    orientation_widget -> SetViewport( 0.0, 0.0, 0.2, 0.2 );
    orientation_widget -> SetEnabled( 1 );
    orientation_widget -> InteractiveOff();


    this -> cell_picker = vtkSmartPointer<PickInteractorStyle>::New();
    this -> cell_picker -> SetDefaultRenderer(this -> renderer);
    this -> cell_picker -> SetCurrentRenderer(this -> renderer);
    this -> cell_picker -> SetMainwindow(this);



    render_window -> GetInteractor() -> SetInteractorStyle(this -> cell_picker);


}



void Mainwindow::init_right_dockwidget(){

 // The lateral dockwidget is filled up
    this -> right_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );

    QWidget * right_dockwidget_container = new QWidget(this);
    QVBoxLayout * right_dockwidget_container_layout = new QVBoxLayout();
    QWidget * selection_widget = new QWidget(this);
    QHBoxLayout * selection_widget_layout = new QHBoxLayout(selection_widget);

    this -> select_facets_button = new QPushButton(selection_widget);
    this -> select_points_button = new QPushButton(selection_widget);
    this -> edit_selection_button = new QPushButton(this);

    selection_widget_layout -> addWidget(select_facets_button);
    selection_widget_layout -> addWidget(select_points_button);

    this -> select_facets_button -> setText("Select facets");
    this -> select_points_button -> setText("Select points");
    this -> edit_selection_button -> setText("Edit selection");

    right_dockwidget_container -> setLayout(right_dockwidget_container_layout);
    
    right_dockwidget_container_layout -> addWidget(this -> prop_table);
    right_dockwidget_container_layout -> addWidget(selection_widget);
    right_dockwidget_container_layout -> addWidget(this -> edit_selection_button);
    right_dockwidget_container_layout -> addWidget(this -> log_console );

    this -> right_dockwidget -> setWidget(right_dockwidget_container);
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> right_dockwidget);

    connect(this -> select_facets_button,SIGNAL(clicked(bool)),this, SLOT(select_facets()));
    connect(this -> select_points_button,SIGNAL(clicked(bool)),this, SLOT(select_points()));
    connect(this -> edit_selection_button,SIGNAL(clicked(bool)),this, SLOT(edit_selection()));

    this -> select_points_button -> setEnabled(1);
    this -> select_facets_button -> setEnabled(0);
    
}

void Mainwindow::select_facets(){

 PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonDown();
 PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonUp();

 this -> select_facets_button -> setEnabled(0);
 this -> select_points_button -> setEnabled(1);

 std::string opening_line = "Switching to facet selection\n";

 std::string closing_line(opening_line.length() - 1, '#');
 closing_line.append("\n");

 this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
 this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
 this -> log_console -> appendPlainText(QString::fromStdString(closing_line));


 this -> facet_selection_mode = true;

};

void Mainwindow::select_points(){

    PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonDown();
    PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonUp();

    this -> select_facets_button -> setEnabled(1);
    this -> select_points_button -> setEnabled(0);

    std::string opening_line = "Switching to vertex selection\n";

    std::string closing_line(opening_line.length() - 1, '#');
    closing_line.append("\n");

    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
    this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
    this -> facet_selection_mode = false;

};



void Mainwindow::set_action_status(bool enabled, QAction * action) {
    action -> setEnabled(enabled);
}

void Mainwindow::update_GUI_changed_prop() {

    // The status bar is updated depending on whether a shape model remains
    if (this -> wrapped_shape_data.size() == 0 ) {
        this -> statusBar() -> showMessage("Ready");
    }

    else if (this -> wrapped_shape_data .size() > 0) {

        int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
        std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

        auto active_shape_polydata  =  this -> wrapped_shape_data[name] -> get_polydata();
        auto N_vertices = active_shape_polydata -> GetPoints() -> GetNumberOfPoints();
        auto N_facets = active_shape_polydata -> GetNumberOfCells();
        std::string message("Facets : " + std::to_string(N_facets) + " Vertices: " + std::to_string(N_vertices));

        this -> statusBar() -> showMessage(QString::fromStdString(message));

    }

    this -> update_actions_availability();

}


void Mainwindow::open_settings_window() {

    SettingsWindow settings_window(this);
    settings_window.exec();

}

void Mainwindow::open_radar_window(){


    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute Simulated Radar Observations", "You must first load a shape model in order to compute radar observations");
        return;
    }
    else{
        RadarWindow radar_window(this);
        radar_window.exec();
    }

}

void Mainwindow::open_lightcurve_window(){

    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute simulated light curve", "You must first load a shape model in order to compute lightcurves");
        return;
    }
    else{
        LCWindow lightcurve_window(this);
        lightcurve_window.exec();
    }
}


void Mainwindow::open_compute_surface_pgm_window(){

    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute Surface PGM", "You must first load a shape model in order to compute/load a surface pgm");
        return;
    }
    else{

        SurfacePGMWindow surface_pgm_window(this);
        surface_pgm_window.exec();
    }

}


void Mainwindow::createActions() {

    this -> open_settings_window_action = new QAction(tr("Preferences"), this);
    this -> open_settings_window_action -> setStatusTip(tr("Open settings window"));
    connect(this -> open_settings_window_action, &QAction::triggered, this, &Mainwindow::open_settings_window);

    this -> save_shape_action = new QAction(tr("Save Shape"), this);
    this -> save_shape_action -> setStatusTip(tr("Save vtkPolyData to OBJ shape model"));
    connect(this -> save_shape_action, &QAction::triggered, this, &Mainwindow::save_shape);

    this -> add_shape_action = new QAction(tr("Load Shape"), this);
    this -> add_shape_action -> setStatusTip(tr("Load obj file holding the facet/vertex description of a shape of interest"));
    connect(this -> add_shape_action, &QAction::triggered, this, &Mainwindow::add_shape);

    this -> clear_console_action = new QAction(tr("Clear Log"), this);
    this -> clear_console_action -> setStatusTip(tr("Clears the log console"));
    connect(this -> clear_console_action, &QAction::triggered, this, &Mainwindow::clear_console);

    this -> open_compute_yorp_window_action = new QAction(tr("Compute YORP Fourier Coefficients"), this);
    this -> open_compute_yorp_window_action -> setStatusTip(tr("Computes the Fourier decomposition of YORP force/torques"));
    connect(this -> open_compute_yorp_window_action, &QAction::triggered, this, &Mainwindow::open_compute_yorp_window);


    this -> open_compute_sharm_window_action = new QAction(tr("Compute Gravity Spherical Harmonics"), this);
    this -> open_compute_sharm_window_action -> setStatusTip(tr("Computes the spherical harmonics coefficients of the exterior gravity field"));
    connect(this -> open_compute_sharm_window_action, &QAction::triggered, this, &Mainwindow::open_compute_sharm_window);


    this -> open_radar_window_action = new QAction(tr("Generate Simulated Radar Observations"), this);
    this -> open_radar_window_action -> setStatusTip(tr("Generates simulated range/range-rate observations emulating a doppler radar"));
    connect(this -> open_radar_window_action, &QAction::triggered, this, &Mainwindow::open_radar_window);


    this -> open_lightcurve_window_action = new QAction(tr("Generate Simulated Light Curve"), this);
    this -> open_lightcurve_window_action -> setStatusTip(tr("Generates simulated light curve"));
    connect(this -> open_lightcurve_window_action, &QAction::triggered, this, &Mainwindow::open_lightcurve_window);


    this -> align_shape_action = new QAction(tr("Align Shape"), this);
    this -> align_shape_action -> setStatusTip(tr("Align selected shape model with barycenter/principal axis"));
    connect(this -> align_shape_action, &QAction::triggered, this, &Mainwindow::align_shape);


    this -> save_console_action = new QAction(tr("Save Log"), this);
    this -> save_console_action -> setStatusTip(tr("Saves log console to a file"));
    connect(this -> save_console_action, &QAction::triggered, this, &Mainwindow::save_console);


    this -> open_mass_properties_window_action = new QAction(tr("Compute Mass Properties"), this);
    this -> open_mass_properties_window_action -> setStatusTip(tr("Compute mass properties of the selected shape"));
    connect(this -> open_mass_properties_window_action, &QAction::triggered, this, &Mainwindow::open_mass_properties_window);


    this -> open_compute_surface_pgm_window_action = new QAction(tr("Compute/Load Surface PGM"), this);
    this -> open_compute_surface_pgm_window_action -> setStatusTip(tr("Evaluate or loads a surface Polyhedron Gravity Model"));
    connect(this -> open_compute_surface_pgm_window_action, &QAction::triggered, this, &Mainwindow::open_compute_surface_pgm_window);

    this -> open_rendering_properties_window_action =new QAction(tr("Rendering Properties"), this);
    this -> open_rendering_properties_window_action -> setStatusTip(tr("Open window enabling one to change the rendering properties"));
    connect(this -> open_rendering_properties_window_action, &QAction::triggered, this, &Mainwindow::open_rendering_properties_window);

    this -> open_select_mapper_window_action = new QAction(tr("Set Results Overlay"), this);
    this -> open_select_mapper_window_action -> setStatusTip(tr("Open window enabling one to set a mapper to represent surface data for a given shape"));
    connect(this -> open_select_mapper_window_action, &QAction::triggered, this, &Mainwindow::open_select_mapper_window);


}



void Mainwindow::open_select_mapper_window(){

    SelectMapperWindow select_mapper_window(this);
    select_mapper_window.exec();
}



void Mainwindow::open_mass_properties_window(){


    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute Mass Properties", "You must first load a shape model in order to evaluate mass properties");
        return;
    }
    else{
        MassPropertiesWindow mass_properties_window(this);
        mass_properties_window.exec();
    }

    
}



void Mainwindow::update_actions_availability() {

    if (this -> wrapped_shape_data.size() == 0){

        this -> align_shape_action -> setEnabled(false);
        this -> save_shape_action -> setEnabled(false);
    }

    else{

        int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
        std::string name = this -> prop_table -> item(selected_row_index, 0)-> text() .toStdString();

        this -> align_shape_action -> setEnabled(true);
        this -> save_shape_action -> setEnabled(true);
        
    }

}



void Mainwindow::open_compute_yorp_window(){


    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute YORP Fourier Coefficients", "You must first load a shape model in order to compute YORP Fourier Coefficients");
        return;
    }
    else{
        YORPWindow yorp_window(this);
        yorp_window.exec();
    }

}


void Mainwindow::open_compute_sharm_window(){
    if(this -> wrapped_shape_data .size() == 0){
        QMessageBox::warning(this, "Compute Gravity Spherical Harmonics", "You must first load a shape model in order to compute Gravity Spherical Harmonics");
        return;
    }

    else{

        SHARMWindow sharm_window(this);
        sharm_window.exec();
    }
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



void Mainwindow::save_shape(){


    std::string default_name;
    int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
    std::string name = this -> prop_table -> item(selected_row_index, 0)-> text() .toStdString();

    default_name = name;
    
    QString fileName = QFileDialog::getSaveFileName(this,tr("Save shape"), QString::fromStdString(default_name), tr("Wavefront file (*.obj)"));

    if (fileName.isEmpty() == false) {
       int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
       std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

       vtkSmartPointer<SBGATObjWriter> writer = vtkSmartPointer<SBGATObjWriter>::New();

       writer -> SetInputData( this -> wrapped_shape_data[name] -> get_polydata());


       writer -> SetFileName(fileName.toStdString().c_str());
       writer -> Update();
       this -> prop_table ->setItem(selected_row_index, 1, new QTableWidgetItem(""));

   }

}

void Mainwindow::add_shape() {

    QString fileName = QFileDialog::getOpenFileName(this,tr("Load shape"), "~/", tr("Wavefront file (*.obj)"));

    if (fileName.isEmpty() == false) {

        QStringList items_length_unit;
        items_length_unit << tr("Meters") << tr("Kilometers");
        bool ok;
        QString length_unit = QInputDialog::getItem(this, tr("Length Unit Of Loaded Shape Model:"),tr("Shape Units:"), items_length_unit, 0, false, &ok);
        
        if (!ok){
            return;
        }

        QMessageBox::StandardButton enforce_centering_aligment = QMessageBox::question(this, "Shape Alignment", "Force centering on barycenter and principal axes alignment?", QMessageBox::No|QMessageBox::Yes);
        double scaling_factor;
        length_unit == "Meters" ?  scaling_factor = 1 : scaling_factor = 1000;

        std::stringstream ss;
        ss.str(std::string());

        std::string opening_line = "### Loading shape ###";
        this -> log_console -> appendPlainText(QString::fromStdString(opening_line));
        this -> log_console -> appendPlainText(QString::fromStdString("- Loading shape from ") + fileName);

        std::chrono::time_point<std::chrono::system_clock> start, end;

        start = std::chrono::system_clock::now();

        // The name of the shape model is extracted from the path
        int dot_index = fileName.lastIndexOf(".");
        int slash_index = fileName.lastIndexOf("/");
        std::string name = (fileName.toStdString()).substr(slash_index + 1 , dot_index - slash_index - 1);
        std::string basic_name = name;

            // A new ModelDataWrapper is created and stored under the name of the shape model
        std::shared_ptr<ModelDataWrapper> model_data = std::make_shared<ModelDataWrapper>();

            // Reading
        vtkNew<vtkOBJReader> reader;
        reader -> SetFileName(fileName.toStdString().c_str());

        reader -> Update();


            // Scaling
        vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
        transform -> Scale(scaling_factor,scaling_factor,scaling_factor);

        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter -> SetInputConnection(reader -> GetOutputPort());
        transformFilter -> SetTransform(transform);
        transformFilter -> SetOutputPointsPrecision ( vtkAlgorithm::DesiredOutputPrecision::DOUBLE_PRECISION );

        transformFilter -> Update();



        vtkSmartPointer<vtkIdFilter> idFilter =
        vtkSmartPointer<vtkIdFilter>::New();
        idFilter -> SetInputConnection(transformFilter->GetOutputPort());
        idFilter -> SetIdsArrayName("OriginalIds");
        idFilter -> Update();

        vtkSmartPointer<vtkDataSetSurfaceFilter> surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
        surfaceFilter -> SetInputConnection(idFilter -> GetOutputPort());
        surfaceFilter -> Update();



        vtkSmartPointer<vtkPolyDataNormals> normals =vtkSmartPointer<vtkPolyDataNormals>::New();
        normals -> SetInputConnection(surfaceFilter-> GetOutputPort());
        normals -> SplittingOff();
        normals -> ConsistencyOn();
        normals -> ComputeCellNormalsOff();
        normals -> ComputePointNormalsOn();
        normals -> Update();

        // Create a PolyData
        vtkSmartPointer<vtkPolyData> polygonPolyData = normals -> GetOutput();

            // The camera is moved to be adjusted to the new shape
        vtkSmartPointer<SBGATMassProperties> center_of_mass_filter =
        vtkSmartPointer<SBGATMassProperties>::New();

        center_of_mass_filter -> SetInputData(polygonPolyData);
        center_of_mass_filter -> Update();
        
        this -> renderer -> GetActiveCamera() -> SetPosition(0, 0, 10 * std::cbrt(3./4. / arma::datum::pi * center_of_mass_filter -> GetVolume()) );

            // Create a mapper and actor
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();

        mapper -> SetInputData(polygonPolyData);

        mapper -> ScalarVisibilityOff();

        vtkSmartPointer<vtkActor> actor =vtkSmartPointer<vtkActor>::New();
        actor -> SetMapper(mapper);

            // Visualize
        this -> renderer -> AddActor(actor);


        // Create the tree
        vtkSmartPointer<vtkKdTreePointLocator> tree =  vtkSmartPointer<vtkKdTreePointLocator>::New();
        tree -> SetDataSet(polygonPolyData);
        tree -> BuildLocator();

        // Store
        model_data -> set_polydata(polygonPolyData);
        model_data -> set_actor(actor);
        model_data -> set_mapper(mapper);
        model_data -> set_scale_factor(scaling_factor);
        model_data -> set_tree(tree);


            // The ModelDataWrapper pointer is stored. 
            // If the name is not already taken, nothing special
        unsigned int count = this -> wrapped_shape_data.count(name);

        if(count == 0){
            this -> wrapped_shape_data[name] = model_data;
        }
        else{

            // otherwise, a suffix is added

            while( this -> wrapped_shape_data.find(name) != this -> wrapped_shape_data.end()){
                std::string suffix = "(" + std::to_string(count) + ")";
                name = basic_name + suffix;
                ++count;
            }

            this -> wrapped_shape_data[name] = model_data;
        }


        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        // The shape table is updated to show the newly loaded shape model
        this -> add_prop_to_table_widget(name);



            // The GUI actions are updated
        this -> update_actions_availability();
        this -> update_GUI_changed_prop();

            // The log console displays the name and content of the loaded shape model
        this -> log_console -> appendPlainText(QString::fromStdString("- Loading completed in ")+ QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

        std::string closing_line(opening_line.length() - 1, '#');
        closing_line.append("\n");
        this -> log_console -> appendPlainText(QString::fromStdString(closing_line));

    // The shape that was just added is already selected
        if (enforce_centering_aligment ==  QMessageBox::Yes){
            this -> align_shape();
        }

    }

    this -> qvtkWidget -> GetRenderWindow() -> Modified();

    this -> qvtkWidget -> GetRenderWindow() -> Render();

}


void Mainwindow::align_shape(){

    int selected_row_index = this  -> prop_table -> selectionModel() -> currentIndex().row();
    std::string name = this  -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    vtkSmartPointer<SBGATMassProperties> center_of_mass_filter =
    vtkSmartPointer<SBGATMassProperties>::New();

    center_of_mass_filter -> SetInputData(this -> wrapped_shape_data[name] -> get_polydata());
    center_of_mass_filter -> Update();

    double center[3];
    center_of_mass_filter -> GetCenterOfMass(center);
    arma::mat::fixed<3,3> principal_axes = center_of_mass_filter -> GetPrincipalAxes();
    auto prv = RBK::dcm_to_prv(principal_axes.t());
    double angle = 180 / arma::datum::pi * arma::norm(prv);
    arma::vec axis = arma::normalise(prv);
    
    vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
    transform -> RotateWXYZ(angle,axis.colptr(0));
    
    transform -> Translate( - center[0],  - center[1],  - center[2]);

    vtkSmartPointer<vtkTransformPolyDataFilter> filter =
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    

    filter -> SetInputData(this -> wrapped_shape_data[name] -> get_polydata());
    filter -> SetTransform(transform);
    filter -> Update();

    vtkSmartPointer<vtkPolyDataMapper> mapper = 
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper -> SetInputConnection(filter -> GetOutputPort());

    mapper -> ScalarVisibilityOff();
    this -> wrapped_shape_data[name] -> set_polydata(filter -> GetOutput());

    this -> wrapped_shape_data[name] -> set_mapper(mapper);
    this -> wrapped_shape_data[name] -> get_actor() -> SetMapper(mapper);
    this -> get_renderer() -> RemoveActor2D(this -> wrapped_shape_data[name] -> get_colorbar_actor());


    // Should update normals and kd tree here
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals -> SetInputData(this -> wrapped_shape_data[name] -> get_polydata());
    normals -> SplittingOff();
    normals -> ConsistencyOn();
    normals -> ComputeCellNormalsOff();
    normals -> ComputePointNormalsOn();
    normals -> Update();

    this -> wrapped_shape_data[name] -> get_polydata() -> GetPointData() -> GetArray("Normals")-> ShallowCopy(normals -> GetOutput() -> GetPointData() -> GetArray("Normals"));

    this -> wrapped_shape_data[name] -> get_polydata() -> GetPointData() -> Modified();
    this -> wrapped_shape_data[name] -> get_polydata() -> Modified();

    this -> qvtkWidget -> GetRenderWindow() -> Render();
    this -> prop_table -> setItem(selected_row_index, 1, new QTableWidgetItem("Modified"));

    std::string displayed_line = "- Translated " + name + " barycenter to (0,0,0) and aligned principal axes with visualization frame";
    
    std::string closing_line(displayed_line.length() - 1, '#');
    closing_line.append("\n");

    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
    this -> log_console -> appendPlainText(QString::fromStdString(displayed_line));
    this -> log_console -> appendPlainText(QString::fromStdString(closing_line));
    

}


void Mainwindow::add_prop_to_table_widget(std::string name) {


    QTableWidgetItem * nameItem = new QTableWidgetItem(QString::fromStdString(name));
    nameItem -> setFlags(Qt::ItemIsSelectable | Qt::ItemIsEnabled);

    this -> prop_table -> insertRow(this -> prop_table -> rowCount());
    this -> prop_table -> setItem(this -> prop_table -> rowCount() - 1, 0, nameItem);

    QWidget * toggle_visibility_button_container = new QWidget(this -> prop_table -> cellWidget(this -> prop_table -> rowCount() - 1, 0));
    QHBoxLayout* toggle_visibility_button_layout = new QHBoxLayout(toggle_visibility_button_container);
    QPushButton * toggle_visibility_button = new QPushButton(toggle_visibility_button_container);
    toggle_visibility_button -> setText("Hide");
    toggle_visibility_button -> setProperty("name", QVariant(QString::fromStdString(name)));
    toggle_visibility_button_layout -> addWidget(toggle_visibility_button);
    toggle_visibility_button_layout -> setAlignment(Qt::AlignCenter);
    toggle_visibility_button_layout -> setContentsMargins(0, 0, 0, 0);
    toggle_visibility_button_container -> setLayout(toggle_visibility_button_layout);
    this -> prop_table -> setCellWidget(this -> prop_table -> rowCount() - 1, 2, toggle_visibility_button_container);

    QWidget * erase_button_container = new QWidget(this -> prop_table -> cellWidget(this -> prop_table -> rowCount() - 1, 0));
    QHBoxLayout* layout = new QHBoxLayout(erase_button_container);
    QPushButton * erase_shape_button = new QPushButton(erase_button_container);
    erase_shape_button -> setText("Delete");
    erase_shape_button -> setProperty("name", QVariant(QString::fromStdString(name)));
    layout -> addWidget(erase_shape_button);
    layout -> setAlignment(Qt::AlignCenter);
    layout -> setContentsMargins(0, 0, 0, 0);
    erase_button_container -> setLayout(layout);
    this -> prop_table -> setCellWidget(this -> prop_table -> rowCount() - 1, 3, erase_button_container);



    // The check box is connected to the proper slot
    connect(toggle_visibility_button, SIGNAL(clicked(bool)), this, SLOT(toggle_prop_visibility()));

    // The check box is connected to the proper slot
    connect(erase_shape_button, SIGNAL(clicked(bool)), this, SLOT(remove_prop()));

    // This prop is selected
    this -> prop_table -> selectRow(this -> prop_table -> rowCount() - 1);

}


void Mainwindow::toggle_prop_visibility() {

 int selected_row_index = this -> prop_table -> selectionModel() -> currentIndex().row();
 std::string name = this -> prop_table -> item(selected_row_index, 0) -> text() .toStdString();

    // Showing/hiding small body shape model actor
     QPushButton * senderObj = qobject_cast<QPushButton*>(sender()); // This will give Sender object

  // This will give obejct name for above it will give "A", "B", "C"
     QString state_string = senderObj -> text(); 

     if (state_string == "Hide") {

        this -> wrapped_shape_data[name] -> get_actor() -> VisibilityOff();

        senderObj -> setText("Show"); 
    }

    else {
        this -> wrapped_shape_data[name] -> get_actor() -> VisibilityOn();

        senderObj -> setText("Hide"); 

    }

    

    // The Render window is updated
    this -> qvtkWidget -> GetRenderWindow() -> Render();

}


void Mainwindow::remove_prop() {

    QPushButton * button = qobject_cast<QPushButton*>(sender());
    std::string name = button -> property("name") . toString().toStdString();



    // The actor of this shape is removed
    this -> renderer -> RemoveActor(this -> wrapped_shape_data[name] -> get_actor());
    this -> renderer -> RemoveActor2D(this -> wrapped_shape_data[name] -> get_colorbar_actor());

    // This emulated left click will lead to the removal of any facet-highlighting actor that could be remaining


    // The data wrapper is removed
    this -> wrapped_shape_data.erase(name);


    // The corresponding row in the table widget is removed
    // This will trigger the corresponding signal/slot mechanism updating the GUI
    for (int i = 0; i < this -> prop_table -> rowCount(); ++i) {
        if (this -> prop_table -> item(i, 0) -> text() == QString::fromStdString(name)) {
            this -> prop_table -> removeRow(i);

            break;
        }
    }


    PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonDown();
    PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> OnLeftButtonUp();


    this -> update_actions_availability();

    // The Render window is updated
    this -> qvtkWidget -> GetRenderWindow() -> Render();

}




void Mainwindow::open_rendering_properties_window(){

    RenderingPropertiesWindow * rendering_properties_window = new RenderingPropertiesWindow(this);

    connect(this,SIGNAL(prop_removed_signal()),rendering_properties_window,SLOT(prop_removed_slot()));
    connect(this,SIGNAL(prop_added_signal()),rendering_properties_window,SLOT(prop_added_slot()));
    
    rendering_properties_window -> show();
}   

void Mainwindow::createMenus() {

    this -> SettingsMenu = this -> menuBar() -> addMenu(tr("&Settings"));
    this -> SettingsMenu -> addAction(this -> open_settings_window_action);

    this -> SmallBodyMenu = this -> menuBar() -> addMenu(tr("&Shape"));
    this -> SmallBodyMenu -> addAction(this -> add_shape_action);
    this -> SmallBodyMenu -> addAction(this -> save_shape_action);
    this -> SmallBodyMenu -> addSeparator();
    this -> SmallBodyMenu -> addAction(this -> align_shape_action);

    

    this -> ObservationsMenu = menuBar() -> addMenu(tr("&Observations"));
    this -> ObservationsMenu -> addAction(this -> open_radar_window_action);
    this -> ObservationsMenu -> addAction(this -> open_lightcurve_window_action);



    this -> AnalysesMenu = menuBar() -> addMenu(tr("&Analyses"));
    this -> AnalysesMenu -> addAction(this -> open_mass_properties_window_action);
    this -> AnalysesMenu -> addAction(this -> open_compute_yorp_window_action);
    this -> AnalysesMenu -> addAction(this -> open_compute_sharm_window_action);
    this -> AnalysesMenu -> addAction(this -> open_compute_surface_pgm_window_action);

    this -> ResultsMenu = menuBar() -> addMenu(tr("&Visualization"));
    this -> ResultsMenu -> addAction(this -> open_select_mapper_window_action);
    this -> ResultsMenu -> addAction(this -> open_rendering_properties_window_action);

    this -> ConsoleMenu = menuBar() -> addMenu(tr("&Console"));
    this -> ConsoleMenu -> addAction(this -> clear_console_action);
    this -> ConsoleMenu -> addAction(this -> save_console_action);

}




DataMap Mainwindow::get_wrapped_shape_data() const{
    return this -> wrapped_shape_data;
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

bool Mainwindow::get_selection_mode()const{
    return this -> facet_selection_mode;
}


void Mainwindow::edit_selection(){

    int selection_size = PickInteractorStyle::SafeDownCast(this -> qvtkWidget -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle()) -> GetSelectionSize();

    if (selection_size != 0){
        if (this -> facet_selection_mode){

            // FacetEditionWindow facet_edition_window(this);
            // facet_edition_window.show();
        }
        else{

            VertexEditionWindow vertex_edition_window(this);
            vertex_edition_window.exec();

        }
    }

}




