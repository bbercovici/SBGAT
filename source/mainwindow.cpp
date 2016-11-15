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

    this -> qvtkWidget -> GetRenderWindow() -> Render();

}

void Mainwindow::setupUi() {
    this -> resize(990, 583);
    this -> lateral_dockwidget = new QDockWidget(this);
    this -> qvtkWidget = new QVTKWidget(this);

    lateral_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );
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
    SelectedPointWidget * selected_point_widget =
        dynamic_cast<SelectedPointWidget *>(this -> lateral_dockwidget -> widget());

    ComputePGMWidget * compute_PGM_widget =
        dynamic_cast<ComputePGMWidget *>(this -> lateral_dockwidget -> widget());

    if (selected_point_widget != nullptr) {
        selected_point_widget -> close();
    }
    else if (compute_PGM_widget != nullptr) {
        compute_PGM_widget -> close();
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

        vtkSmartPointer<vtkExtractEdges> extractEdges =
            vtkSmartPointer<vtkExtractEdges>::New();
        extractEdges -> SetInputData(read_polydata);
        extractEdges -> Update();


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
        this -> set_action_status(true, shapeColorAct);
        this -> set_action_status(true, vertexVisibilityAct);
        this -> set_action_status(true, selectPointAct);
        this -> set_action_status(true, openComputePGMWidgetAct);
        this -> set_action_status(true, saveAct);

        // The topology of the shape is constructed
        read_polydata -> BuildLinks();

        // The interactor is set up and connected to the shape model being displayed
        vtkSmartPointer<InteractorStyle> style =
            vtkSmartPointer<InteractorStyle>::New();

        this -> renderWindowInteractor -> SetInteractorStyle( style );

        style -> set_mainwindow(this);
        style -> set_all_points_polydata(read_polydata);

        this -> qvtkWidget -> GetRenderWindow() -> Render();

    }

}


void Mainwindow::open() {

    // The file name is retrieved from the output of the QFileDialog window
    QString fileName = QFileDialog::getOpenFileName(this,
                       tr("Open File"), "", tr("OBJ jile ( *.obj)"));


    if (!fileName.isEmpty()) {
        this -> close_lateral_dockwidget();

        // A .obj reader is created, connected to the file and fetched into a vtkPolyData structure
        vtkSmartPointer<vtkOBJReader> reader =
            vtkSmartPointer<vtkOBJReader>::New();
        reader -> SetFileName(fileName.toStdString().c_str());
        reader -> Update();
        vtkSmartPointer<vtkPolyData> read_polydata_without_id = reader -> GetOutput();
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

}

void Mainwindow::save() {

    // The save path is queried
    QString fileName = QFileDialog::getSaveFileName(this,
                       tr("Save File"));
    if (!fileName.isEmpty()) {

        // the shape model currently displayed in the Rendering window is saved to a .obj file
        // by means of the appropriate writer class
        vtkSmartPointer<vtkOBJExporter> exporter =
            vtkSmartPointer<vtkOBJExporter>::New();
        exporter -> SetFilePrefix(fileName.toStdString().c_str());
        exporter -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow() );
        exporter -> Update();
        exporter -> Write();

    }
}


void Mainwindow::createActions() {
    openAct = new QAction(tr("&Open"), this);
    openAct -> setShortcuts(QKeySequence::Open);
    openAct -> setStatusTip(tr("Open an existing .vtp file"));
    connect(openAct, &QAction::triggered, this, &Mainwindow::open);


    saveAct = new QAction(tr("&Save"), this);
    saveAct -> setShortcuts(QKeySequence::Save);
    saveAct -> setStatusTip(tr("Save the current vtkPolyDataObject to a .vtp file"));
    connect(saveAct, &QAction::triggered, this, &Mainwindow::save);
    saveAct -> setDisabled(true);

    selectPointAct = new QAction(tr("Select Points"), this);
    selectPointAct -> setStatusTip(tr("Select a group of vertices and operates on them"));
    connect(selectPointAct, &QAction::triggered, this, &Mainwindow::select);
    selectPointAct -> setDisabled(true);

    openComputePGMWidgetAct = new QAction(tr("Compute Polyhedron Gravity Model"), this);
    openComputePGMWidgetAct -> setStatusTip(tr("Open widget enabling computation of the Polyhedron Gravity Model of the displayed shape assuming a constant density distribution"));
    connect(openComputePGMWidgetAct, &QAction::triggered, this, &Mainwindow::open_compute_pgm_widget);
    openComputePGMWidgetAct -> setDisabled(true);



    shapeColorAct = new QAction(tr("Shape Color"), this);
    shapeColorAct -> setStatusTip(tr("Set the shape actor color"));
    connect(shapeColorAct, &QAction::triggered, this, &Mainwindow::set_shape_color);
    shapeColorAct -> setDisabled(true);


    backgroundColorAct = new QAction(tr("Background Color"), this);
    backgroundColorAct -> setStatusTip(tr("Set the background color"));
    connect(backgroundColorAct, &QAction::triggered, this, &Mainwindow::set_background_color);


    resetAct = new QAction(tr("Clear all"), this);
    resetAct -> setStatusTip(tr("Clear all windows and loaded shape models"));
    connect(resetAct, &QAction::triggered, this, &Mainwindow::clear_all);


    vertexVisibilityAct = new QAction(tr("Show Vertices"), this);
    vertexVisibilityAct -> setStatusTip(tr("Display or Hide Vertices"));
    vertexVisibilityAct -> setCheckable(true);
    vertexVisibilityAct -> setChecked(true);
    vertexVisibilityAct -> setDisabled(true);
    connect(vertexVisibilityAct, &QAction::triggered, this, &Mainwindow::change_vertex_visibility);


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


void Mainwindow::change_vertex_visibility() {
    // The first actor is shown/hiddenr
    vtkSmartPointer<vtkActor> point_cloud_actor = vtkActor::SafeDownCast(
                this -> renderer -> GetActors() -> GetItemAsObject(0));

    // The visibility state is retrieved and inverted
    int visibility = point_cloud_actor -> GetVisibility();

    point_cloud_actor -> SetVisibility(!visibility);
    point_cloud_actor -> Modified();
    this -> qvtkWidget -> GetRenderWindow() -> Render();
}

void Mainwindow::createMenus() {
    fileMenu = menuBar()->addMenu(tr("&File"));

    fileMenu -> addAction(openAct);
    fileMenu -> addAction(saveAct);
    fileMenu -> addSeparator();
    fileMenu -> addAction(resetAct);


    OperationMenu = menuBar() -> addMenu(tr("&Operation"));
    OperationMenu -> addAction(selectPointAct);
    OperationMenu -> addAction(openComputePGMWidgetAct);

    ViewMenu = menuBar() -> addMenu(tr("&Shape Graphic Properties"));
    ViewMenu -> addAction(shapeColorAct);
    ViewMenu -> addAction(backgroundColorAct);
    ViewMenu -> addSeparator();
    ViewMenu -> addAction(vertexVisibilityAct);

}

vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}

void Mainwindow::clear_all() {
    this -> remove_actors();
    this -> close_lateral_dockwidget();
    this -> qvtkWidget -> GetRenderWindow() -> Render();


}

void Mainwindow::open_compute_pgm_widget() {

    ComputePGMWidget * compute_PGM_widget = new ComputePGMWidget(this);

    this ->  lateral_dockwidget -> setWidget(compute_PGM_widget);
    this ->  lateral_dockwidget -> show();



}

