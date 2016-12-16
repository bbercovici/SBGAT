#include "Mainwindow.hpp"

// shortcut to interactor modes
#define INTERACTOR_IS_ORIENT 0
#define INTERACTOR_IS_SELECT 1

Mainwindow::Mainwindow() {
    this -> setupUi();

    // A VTK renderer is created and linked with the qvtk widget
    this -> renderer = vtkSmartPointer<vtkRenderer>::New();
    this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(this -> renderer);

    this -> renderer -> SetGradientBackground (true);
    this -> renderer -> SetBackground (0.5, 0.5, 1);

    vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
    this -> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    this -> renderWindowInteractor -> SetPicker(areaPicker);
    this -> renderWindowInteractor -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow());

    vtkSmartPointer<vtkAxesActor> axes =
        vtkSmartPointer<vtkAxesActor>::New();

    this -> widget =
        vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    widget -> SetOrientationMarker( axes );
    widget -> SetInteractor( this -> renderWindowInteractor );
    widget -> SetViewport( 0.0, 0.0, 0.2, 0.2 );
    widget -> SetEnabled( 1 );
    widget -> InteractiveOff();

    // Setting this -> asteroid to nullptr here prevents the
    // program from crashing when delete is called on this -> asteroid
    // before the first asteroid is loaded
    this -> asteroid = nullptr;

    this -> qvtkWidget -> GetRenderWindow() -> Render();
}

void Mainwindow::setupUi() {
    this -> resize(990, 583);
    this -> lateral_dockwidget = new QDockWidget(this);
    this -> qvtkWidget = new QVTKWidget(this);
    this -> status_bar = new QStatusBar(this);

    this -> setStatusBar(this -> status_bar);
    this -> statusBar() -> showMessage("Ready");

    this -> lateral_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );
    this -> lateral_dockwidget -> hide();

    createActions();
    createMenus();

    this -> setCentralWidget(qvtkWidget);
    this -> setWindowTitle(QStringLiteral("SBGAT (WIP)"));
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> lateral_dockwidget);

    disableGLHiDPI(this -> qvtkWidget -> winId());

    this -> show();

}

void Mainwindow::close_lateral_dockwidget() {
    ModifyAreaWidget * modify_area_widget =
        dynamic_cast<ModifyAreaWidget *>(this -> lateral_dockwidget -> widget());

    ComputePGMWidget * compute_PGM_widget =
        dynamic_cast<ComputePGMWidget *>(this -> lateral_dockwidget -> widget());

    ShapeInfoWidget * shape_info_widget =
        dynamic_cast<ShapeInfoWidget *>(this -> lateral_dockwidget -> widget());


    SetInputScalingWidget * set_input_scaling_widget =
        dynamic_cast<SetInputScalingWidget *>(this -> lateral_dockwidget -> widget());



    if (modify_area_widget != nullptr) {
        modify_area_widget -> close();
    }

    else if (compute_PGM_widget != nullptr) {
        compute_PGM_widget -> close();
    }

    else if (shape_info_widget != nullptr) {
        shape_info_widget -> close();
    }

    else if (set_input_scaling_widget != nullptr) {
        set_input_scaling_widget -> close();
    }


}

vtkSmartPointer<vtkRenderWindowInteractor> Mainwindow::get_render_window_interactor() {
    return this -> renderWindowInteractor;
}

void Mainwindow::set_action_status(bool enabled, QAction * action) {
    action -> setEnabled(enabled);
}

void Mainwindow::load_obj(vtkSmartPointer<vtkPolyData> read_polydata_without_id) {

    if (read_polydata_without_id != NULL) {

        // The polydata is fed to an IDFilter.
        vtkSmartPointer<vtkIdFilter>  id_filter = vtkSmartPointer<vtkIdFilter>::New();
        id_filter -> SetIdsArrayName("ids");
        id_filter -> SetInputData(read_polydata_without_id);
        id_filter -> PointIdsOn();
        id_filter -> Update();

        vtkSmartPointer<vtkPolyData> read_polydata = id_filter -> GetPolyDataOutput();

        // Create a mapper and actor to represent the shape model
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper -> SetInputData(read_polydata);
        vtkSmartPointer<vtkActor> shape_actor = vtkSmartPointer<vtkActor>::New();
        shape_actor -> SetMapper(mapper);
        shape_actor -> GetMapper() -> ScalarVisibilityOff();
        shape_actor -> GetProperty() -> SetColor(1, 1, 1);
        shape_actor -> GetProperty() -> SetPointSize(3);


        // Remove any already existing actor from the rendering window
        this -> remove_actors();
        actor_vector.push_back(shape_actor);

        // The new actor is added to the scene
        this -> renderer -> AddActor(shape_actor);

        // The actions that were until now disabled are enabled
        this -> set_action_status(true, this -> set_shape_color_action);
        this -> set_action_status(true, this -> modify_shape_action);
        this -> set_action_status(true, this -> open_ComputePGMWidget_action);
        this -> set_action_status(true, this -> open_ShapeInfoWidget_action);
        this -> set_action_status(true, this -> save_action);

        // The topology of the shape is constructed
        read_polydata -> BuildLinks();

        // The interactor is set up and connected to the shape model being displayed
        vtkSmartPointer<InteractorStyle> style =
            vtkSmartPointer<InteractorStyle>::New();

        this -> renderWindowInteractor -> SetInteractorStyle( style );
        style -> set_mainwindow(this);

        // Any previously loaded asteroid is deleted 
        // and replaced by a new one of shape model defined by 
        // the polydata just read in
        delete(this -> asteroid);
        this -> asteroid = new Asteroid(read_polydata, 1000);

        // The camera position is adjusted
        this -> renderer -> GetActiveCamera () -> SetPosition (
            read_polydata -> GetLength() * 1.,
            read_polydata -> GetLength() * 1.,
            read_polydata -> GetLength() * 1.);
        this -> statusBar() -> showMessage("Ready");

        this -> renderer -> ResetCamera();
        this -> qvtkWidget -> GetRenderWindow() -> Render();

    }

}

std::vector<vtkSmartPointer<vtkActor> > Mainwindow::get_actor_vector() {
    return this -> actor_vector;
}

Asteroid * Mainwindow::get_asteroid() {
    return this -> asteroid;
}

void Mainwindow::open() {

    // The file name is retrieved from the output of the QFileDialog window
    QString fileName = QFileDialog::getOpenFileName(this,
                       tr("Open File"), "../resources/", tr("OBJ file ( *.obj)"));

    if (!fileName.isEmpty()) {
        this -> statusBar() -> showMessage("Opening .obj");

        this -> close_lateral_dockwidget();

        // The user indicates which units are used in the read .obj file
        // This widget sets the scaling_factor member of mainwindow
        // to the input value
        SetInputScalingWidget * set_input_scaling_widget = new SetInputScalingWidget(this);
        this ->  lateral_dockwidget -> setWidget(set_input_scaling_widget);
        this ->  lateral_dockwidget -> show();
        set_input_scaling_widget -> exec();

        // A .obj reader is created, connected to the file and fetched into a vtkPolyData structure
        vtkSmartPointer<vtkOBJReader> reader =
            vtkSmartPointer<vtkOBJReader>::New();
        reader -> SetFileName(fileName.toStdString().c_str());
        reader -> Update();

        // A VTKTransformFilter is created and provided
        // with an appropriate scaling transform
        vtkSmartPointer<vtkTransform> transform  =
            vtkSmartPointer<vtkTransform>::New();
        transform -> Scale(this -> scaling_factor, 
            this -> scaling_factor, 
            this -> scaling_factor);

        // A scaling transform is applied to the input polydata so as to have its vertex
        // coordinates expressed in meters
        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
            vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter -> SetInputConnection(reader -> GetOutputPort());
        transformFilter -> SetTransform(transform);
        transformFilter -> Update();
        vtkSmartPointer<vtkPolyData> read_polydata_without_id = transformFilter -> GetOutput();
        
        // The content of the obj file is loaded into SBGAT 
        this -> load_obj(read_polydata_without_id);

    }


}

void Mainwindow::set_background_color() {
    QColorDialog *  palette = new QColorDialog(this);

    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> SetBackground (float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}

void Mainwindow::set_shape_color() {
    QColorDialog *  palette = new QColorDialog(this);

    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> GetActors () -> GetLastActor()
        -> GetProperty() -> SetColor(float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}


void Mainwindow::select() {
    // Allows the interactor to grab props by
    // setting its style mode to selection
    // this is equivalent to pressing the "r" key
    // when the interactor is in ORIENT mode
    InteractorStyle * style = static_cast< InteractorStyle * >(this -> qvtkWidget
                              -> GetRenderWindow() -> GetInteractor() -> GetInteractorStyle());

    style -> set_current_mode(INTERACTOR_IS_SELECT);
    this -> close_lateral_dockwidget();


}

void Mainwindow::save() {

    // The save path is queried
    QString fileName = QFileDialog::getSaveFileName(this,
                       "Save File", "../saved_pgm/", tr("OBJ File(*.obj)"));
    if (!fileName.isEmpty()) {


        // A VTKTransformFilter is created and provided
        // with an appropriate scaling transform
        vtkSmartPointer<vtkTransform> transform  =
            vtkSmartPointer<vtkTransform>::New();
        transform -> Scale(1. / this -> scaling_factor,
                           1. / this -> scaling_factor,
                           1. / this -> scaling_factor);

        vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter =
            vtkSmartPointer<vtkTransformPolyDataFilter>::New();
        transformFilter -> SetInputData(this -> asteroid -> get_polydata());
        transformFilter -> SetTransform(transform);
        transformFilter -> Update();

        vtkSmartPointer<vtkOBJWriter> writer =
            vtkSmartPointer<vtkOBJWriter>::New();
        writer -> SetInputData(transformFilter -> GetOutput() );
        writer -> SetFileName(fileName.toStdString().c_str());
        writer -> Update();
    }
}


void Mainwindow::createActions() {
    this -> load_action = new QAction(tr("&Load"), this);
    this -> load_action -> setShortcuts(QKeySequence::Open);
    this -> load_action -> setStatusTip(tr("Loads an existing .obj shape model"));
    connect(this -> load_action, &QAction::triggered, this, &Mainwindow::open);


    this -> save_action = new QAction(tr("&Save"), this);
    this -> save_action -> setShortcuts(QKeySequence::Save);
    this -> save_action -> setStatusTip(tr("Save the current vtkPolyDataObject to a .vtp file"));
    connect(this -> save_action, &QAction::triggered, this, &Mainwindow::save);
    this -> save_action -> setDisabled(true);

    this -> modify_shape_action = new QAction(tr("Modify Area"), this);
    this -> modify_shape_action -> setStatusTip(tr("Select a group of vertices and provides shape modifying options"));
    connect(this -> modify_shape_action, &QAction::triggered, this, &Mainwindow::select);
    this -> modify_shape_action -> setDisabled(true);

    this -> open_ComputePGMWidget_action = new QAction(tr("Polyhedron Gravity Model"), this);
    this -> open_ComputePGMWidget_action -> setStatusTip(tr("Open widget enabling computation of the Polyhedron Gravity Model of the displayed shape assuming a constant density distribution"));
    connect(this -> open_ComputePGMWidget_action, &QAction::triggered, this, &Mainwindow::open_compute_pgm_widget);
    this -> open_ComputePGMWidget_action -> setDisabled(true);


    this -> open_ShapeInfoWidget_action = new QAction(tr("Show Shape Model Info"), this);
    this -> open_ShapeInfoWidget_action -> setStatusTip(tr("Open widget showing shape model information"));
    connect(this -> open_ShapeInfoWidget_action, &QAction::triggered, this, &Mainwindow::open_shape_info_widget);
    this -> open_ShapeInfoWidget_action -> setDisabled(true);

    this -> set_shape_color_action = new QAction(tr("Shape Color"), this);
    this -> set_shape_color_action -> setStatusTip(tr("Set the shape actor color"));
    connect(this -> set_shape_color_action, &QAction::triggered, this, &Mainwindow::set_shape_color);
    this -> set_shape_color_action -> setDisabled(true);


    this -> set_background_color_action = new QAction(tr("Background Color"), this);
    this -> set_background_color_action -> setStatusTip(tr("Set the background color"));
    connect(this -> set_background_color_action, &QAction::triggered, this, &Mainwindow::set_background_color);


    this -> clear_all_action = new QAction(tr("Clear all"), this);
    this -> clear_all_action -> setStatusTip(tr("Clear all windows and loaded shape models"));
    connect(this -> clear_all_action, &QAction::triggered, this, &Mainwindow::clear_all);

}

void Mainwindow::set_actors_visibility(bool visibility) {
    for (std::vector<vtkSmartPointer<vtkActor > >::iterator iter = actor_vector.begin();
            iter != actor_vector.end(); ++ iter) {
        (*iter) -> SetVisibility(visibility);
    }
}

void Mainwindow::remove_actors() {
    for (std::vector<vtkSmartPointer<vtkActor> >::iterator iter = this -> actor_vector.begin();
            iter != this -> actor_vector.end(); ++iter) {
        this -> get_renderer() -> RemoveActor(*iter);
    }
    actor_vector.clear();
}


void Mainwindow::createMenus() {
    this -> fileMenu = menuBar()->addMenu(tr("&File"));

    this -> fileMenu -> addAction(this -> load_action);
    this -> fileMenu -> addAction(this -> save_action);
    this -> fileMenu -> addSeparator();
    this -> fileMenu -> addAction(this -> clear_all_action);

    this -> ShapeModelMenu = menuBar() -> addMenu(tr("&Shape Model"));
    this -> ShapeModelMenu -> addAction(this -> modify_shape_action);
    this -> ShapeModelMenu -> addAction(this -> open_ComputePGMWidget_action);
    this -> ShapeModelMenu -> addSeparator();
    this -> ShapeModelMenu -> addAction(this -> open_ShapeInfoWidget_action);

    this -> ViewMenu = menuBar() -> addMenu(tr("&Shape Graphic Properties"));
    this -> ViewMenu -> addAction(this -> set_shape_color_action);
    this -> ViewMenu -> addAction(this -> set_background_color_action);
    this -> ViewMenu -> addSeparator();


}

vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}

void Mainwindow::clear_all() {
    // All the actors currently displayed are removed
    this -> remove_actors();

    // The lateral dockwidget that was opened (if any) is closed 
    // and its close() method called
    this -> close_lateral_dockwidget();

    // The current asteroid is destroyed 
    delete(this -> asteroid);
    this -> asteroid = nullptr;

    // The render window is updated
    this -> qvtkWidget -> GetRenderWindow() -> Render();

    // The actions that were until now enabled are disabled
    this -> set_action_status(false, this -> set_shape_color_action);
    this -> set_action_status(false, this -> modify_shape_action);
    this -> set_action_status(false, this -> open_ComputePGMWidget_action);
    this -> set_action_status(false, this -> open_ShapeInfoWidget_action);
    this -> set_action_status(false, this -> save_action);



}

void Mainwindow::open_compute_pgm_widget() {

    this -> close_lateral_dockwidget();
    ComputePGMWidget * compute_PGM_widget = new ComputePGMWidget(this);

    this ->  lateral_dockwidget -> setWidget(compute_PGM_widget);
    this ->  lateral_dockwidget -> show();

}


void Mainwindow::open_shape_info_widget() {

    this -> close_lateral_dockwidget();
    ShapeInfoWidget * shape_info_widget = new ShapeInfoWidget(this);

    this ->  lateral_dockwidget -> setWidget(shape_info_widget);
    this ->  lateral_dockwidget -> show();

}

void Mainwindow::set_scaling_factor(double scaling_factor) {
    this -> scaling_factor = scaling_factor;
}


