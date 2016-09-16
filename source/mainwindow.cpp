#include "mainwindow.h"
#include "utilities.h"

// shortcut to interactor modes
#define INTERACTOR_IS_ORIENT 0
#define INTERACTOR_IS_SELECT 1

MainWindow::MainWindow() {
    this -> setupUi();

    // A VTK renderer is created and linked with the qvtk widget
    this -> renderer = vtkSmartPointer<vtkRenderer>::New();
    this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(this -> renderer);

    this -> renderer -> SetGradientBackground (true);
    this -> renderer -> SetBackground (0.5, 0.5, 1);

    this -> selection_widget_is_open = new bool();
    *this -> selection_widget_is_open = false;


    vtkSmartPointer<vtkAreaPicker> areaPicker = vtkSmartPointer<vtkAreaPicker>::New();
    this -> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    this -> renderWindowInteractor -> SetPicker(areaPicker);
    this -> renderWindowInteractor -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow());

    this -> style =
        vtkSmartPointer<InteractorStyle>::New();
    this -> renderWindowInteractor -> SetInteractorStyle( style );

    this -> style -> set_mainwindow(this);

    this -> qvtkWidget -> GetRenderWindow() -> Render();

    leak_tracker = vtkDebugLeaks::New();


}

void MainWindow::setupUi() {
    this -> resize(990, 583);
    selected_point_dockwidget = new QDockWidget(this);
    palette = new QColorDialog(this);
    qvtkWidget = new QVTKWidget(this);
    pc_editing_widget = new SelectedPointWidget(this);

    selected_point_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );

    // hides the dock widget
    pc_editing_widget -> QDialog::reject();

    createActions();
    createMenus();

    this -> setCentralWidget(qvtkWidget);
    this -> setWindowTitle(QStringLiteral("PDART (WIP)"));
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> selected_point_dockwidget);
    this -> selected_point_dockwidget -> setWidget(pc_editing_widget);

    pc_editing_widget -> hide();

    this -> show();

}

void MainWindow::load_pc(vtkSmartPointer<vtkPolyData> read_polydata_without_id) {
    if (read_polydata_without_id != NULL) {

        // The polydata is fed to an IDFilter.
        vtkSmartPointer<vtkIdFilter>  id_filter = vtkSmartPointer<vtkIdFilter>::New();
        id_filter -> SetIdsArrayName("ids");
        id_filter -> SetInputData(read_polydata_without_id);
        id_filter -> PointIdsOn();
        id_filter -> Update();

        vtkSmartPointer<vtkPolyData> read_polydata = id_filter -> GetPolyDataOutput();

        // Create a mapper and actor to represent the point cloud
        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper -> SetInputData(read_polydata);

        vtkSmartPointer<vtkActor> point_cloud_actor = vtkSmartPointer<vtkActor>::New();

        point_cloud_actor -> SetMapper(mapper);
        point_cloud_actor -> GetMapper()->ScalarVisibilityOff();
        point_cloud_actor -> GetProperty()->SetColor(1, 1, 1);
        point_cloud_actor -> GetProperty()->SetPointSize(3);

        // A convex hull is generated
        vtkSmartPointer<vtkDelaunay3D> delaunay3D =
            vtkSmartPointer<vtkDelaunay3D>::New();
        delaunay3D -> SetInputDataObject (read_polydata);

        // The convex hull is mapped, fed into an actor and displayed
        vtkSmartPointer<vtkDataSetMapper> delaunayMapper =
            vtkSmartPointer<vtkDataSetMapper>::New();
        delaunayMapper -> ScalarVisibilityOff();
        delaunayMapper -> SetInputConnection(delaunay3D -> GetOutputPort());

        vtkSmartPointer<vtkActor> delaunay_actor =
            vtkSmartPointer<vtkActor>::New();
        delaunay_actor -> SetMapper(delaunayMapper);
        delaunay_actor -> GetProperty() -> SetColor(0, 1, 0);

        // Remove any already existing actor from the rendering window
        for (std::vector<vtkSmartPointer<vtkActor>>::iterator iter = this -> actor_vector.begin();
                iter != this -> actor_vector.end(); ++iter) {
            this -> renderer -> RemoveActor(*iter);
        }
        actor_vector.clear();

        // The new actors are added to the scene
        this -> renderer -> AddActor(point_cloud_actor);
        this -> renderer -> AddActor(delaunay_actor);

        // The actions that were until now disabled are enabled
        this -> set_action_status(true, "&Shape Graphic Properties", "Shape Color");
        this -> set_action_status(true, "&Shape Graphic Properties", "Show Vertices");
        this -> set_action_status(true, "&Operation", "Select Points");
        this -> set_action_status(true, "&File", "&Save");


        // The interactor is set up and connected to the shape model being displayed

        style -> set_all_points_polydata(read_polydata);


        // The pointers to both actors are stored
        actor_vector.push_back(point_cloud_actor);
        actor_vector.push_back(delaunay_actor);





    }

}



void MainWindow::load_obj(vtkSmartPointer<vtkPolyData> read_polydata_without_id) {

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
        this -> set_action_status(true, "&Shape Graphic Properties", "Shape Color");
        this -> set_action_status(true, "&Shape Graphic Properties", "Show Vertices");
        this -> set_action_status(true, "&Operation", "Select Points");
        this -> set_action_status(true, "&File", "&Save");

        // The topology of the shape is constructed
        read_polydata -> BuildLinks();

        // The interactor is set up and connected to the shape model being displayed
        style -> set_all_points_polydata(read_polydata);
        

        this -> qvtkWidget -> GetRenderWindow() -> Render();

    }

}

void MainWindow::set_action_status(bool enabled, const std::string & menu_name, const std::string & action_name) {
    QList<QMenu*> menus = this -> menuBar() -> findChildren<QMenu*>();
    for (int i = 0; i < menus.size(); ++i) {
        if ( menus.at(i) -> title().toStdString().compare(menu_name) == 0 ) {
            QList<QAction*> actions = menus.at(i)->actions();
            for (int j = 0; j < actions.size(); ++j) {
                if (actions.at(j) -> text().toStdString().compare(action_name) == 0) {
                    actions.at(j) -> setEnabled(enabled);
                }
            }

        }
    }
}

void MainWindow::open() {

    // The file name is retrieved from the output of the QFileDialog window
    QString fileName = QFileDialog::getOpenFileName(this,
                       tr("Open File"), "", tr("VTK PolyData file or OBJ jile (*.vtp , *.obj)"));


    if (!fileName.isEmpty()) {

        // The file extension is extracted
        std::string filename_std_string = fileName.toStdString();
        std::string file_ex = "";
        size_t i = filename_std_string.rfind('.', filename_std_string.length());

        // if an extension can be extracted
        if (i != std::string::npos) {

            // the extension is first stored in a string
            file_ex = filename_std_string.substr(i + 1, filename_std_string.length() - i);

            // for convenience purposes, the string is associated with an integer value
            // representing the corresponding file type
            int file_type = -1;

            if (file_ex.compare("vtp") == 0) {
                file_type = 0;
            }
            if (file_ex.compare("obj") == 0) {
                file_type = 1;
            }

            // If the file has the right extension, the current job is terminated
            if (file_type > 0) {
                this -> terminate_current_job ();
            }

            switch (file_type) {

            case 0: {
                // A .vtp reader is created, connected to the file and fetched into a vtkPolyData structure
                vtkSmartPointer<vtkXMLPolyDataReader> reader =
                    vtkSmartPointer<vtkXMLPolyDataReader>::New();
                reader -> SetFileName(filename_std_string.c_str());
                reader -> Update();
                vtkSmartPointer<vtkPolyData> read_polydata_without_id = reader -> GetOutput();
                this -> load_pc(read_polydata_without_id);

                break;
            }

            case 1: {

                // A .obj reader is created, connected to the file and fetched into a vtkPolyData structure
                vtkSmartPointer<vtkOBJReader> reader =
                    vtkSmartPointer<vtkOBJReader>::New();
                reader -> SetFileName(filename_std_string.c_str());
                reader -> Update();
                vtkSmartPointer<vtkPolyData> read_polydata_without_id = reader -> GetOutput();
                this -> load_obj(read_polydata_without_id);

                break;
            }

            default:
                break;

            }


        }
    }


}

void MainWindow::set_background_color() {
    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> SetBackground (float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}

void MainWindow::set_shape_color() {
    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> GetActors () -> GetLastActor()
        -> GetProperty() -> SetColor(float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}


void MainWindow::select() {
    // Allows the interactor to grab props by
    // setting its style mode to selection
    // this is equivalent to pressing the "r" key
    // when the interactor is in ORIENT mode
    this -> style -> set_current_mode(INTERACTOR_IS_SELECT);

}

void MainWindow::save() {

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


void MainWindow::createActions() {
    openAct = new QAction(tr("&Open"), this);
    openAct -> setShortcuts(QKeySequence::Open);
    openAct -> setStatusTip(tr("Open an existing .vtp file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::open);


    saveAct = new QAction(tr("&Save"), this);
    saveAct -> setShortcuts(QKeySequence::Save);
    saveAct -> setStatusTip(tr("Save the current vtkPolyDataObject to a .vtp file"));
    connect(saveAct, &QAction::triggered, this, &MainWindow::save);
    saveAct -> setDisabled(true);

    selectPointAct = new QAction(tr("Select Points"), this);
    selectPointAct -> setStatusTip(tr("Select a group of vertices and operates on them"));
    connect(selectPointAct, &QAction::triggered, this, &MainWindow::select);
    selectPointAct -> setDisabled(true);


    shapeColorAct = new QAction(tr("Shape Color"), this);
    shapeColorAct -> setStatusTip(tr("Set the shape actor color"));
    connect(shapeColorAct, &QAction::triggered, this, &MainWindow::set_shape_color);
    shapeColorAct -> setDisabled(true);


    backgroundColorAct = new QAction(tr("Background Color"), this);
    backgroundColorAct -> setStatusTip(tr("Set the background color"));
    connect(backgroundColorAct, &QAction::triggered, this, &MainWindow::set_background_color);


    resetAct = new QAction(tr("Reset"), this);
    resetAct -> setStatusTip(tr("Reset PDART to its start-up state"));
    connect(resetAct, &QAction::triggered, this, &MainWindow::reset);


    vertexVisibilityAct = new QAction(tr("Show Vertices"), this);
    vertexVisibilityAct -> setStatusTip(tr("Display or Hide Vertices"));
    vertexVisibilityAct -> setCheckable(true);
    vertexVisibilityAct -> setChecked(true);
    vertexVisibilityAct -> setDisabled(true);
    connect(vertexVisibilityAct, &QAction::triggered, this, &MainWindow::change_vertex_visibility);
}

void MainWindow::set_actors_visibility(bool visibility) {
    for (std::vector<vtkSmartPointer<vtkActor > >::iterator iter = actor_vector.begin();
            iter != actor_vector.end(); ++ iter) {
        (*iter) -> SetVisibility(visibility);
    }
}


void MainWindow::remove_actors() {
    for (std::vector<vtkSmartPointer<vtkActor> >::iterator iter = this -> actor_vector.begin();
            iter != this -> actor_vector.end(); ++iter) {
        this -> get_renderer() -> RemoveActor(*iter);
    }
    actor_vector.clear();
}


void MainWindow::change_vertex_visibility() {
    // The first actor is shown/hiddenr
    vtkSmartPointer<vtkActor> point_cloud_actor = vtkActor::SafeDownCast(
                this -> renderer -> GetActors() -> GetItemAsObject(0));

    // The visibility state is retrieved and inverted
    int visibility = point_cloud_actor -> GetVisibility();

    point_cloud_actor -> SetVisibility(!visibility);
    point_cloud_actor -> Modified();
    this -> qvtkWidget -> GetRenderWindow() -> Render();
}

void MainWindow::createMenus() {
    fileMenu = menuBar()->addMenu(tr("&File"));

    fileMenu -> addAction(openAct);
    fileMenu -> addAction(saveAct);
    fileMenu -> addSeparator();
    fileMenu -> addAction(resetAct);


    OperationMenu = menuBar() -> addMenu(tr("&Operation"));
    OperationMenu -> addAction(selectPointAct);

    ViewMenu = menuBar() -> addMenu(tr("&Shape Graphic Properties"));
    ViewMenu -> addAction(shapeColorAct);
    ViewMenu -> addAction(backgroundColorAct);
    ViewMenu -> addSeparator();
    ViewMenu -> addAction(vertexVisibilityAct);

}

void MainWindow::terminate_current_job() {
    if (*this -> selection_widget_is_open) {
        // pc_editing_widget -> close() automatically calls pc_editing_widget -> reject();
        this -> pc_editing_widget -> close();
    }
}

vtkSmartPointer<vtkRenderer> MainWindow::get_renderer() {
    return this -> renderer;
}


void MainWindow::reset() {
    this -> terminate_current_job();
    this -> remove_actors();
    this -> pc_editing_widget -> reset();
    this -> style -> reset();
    this -> style -> Delete();



//    the tree methods below do not appear to clear memory
//    this -> pc_editing_widget -> reset();
//    this -> style -> reset();


    this -> qvtkWidget -> GetRenderWindow() -> Render();


}

MainWindow::~MainWindow() {
    delete selection_widget_is_open;
}