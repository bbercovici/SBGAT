#include "mainwindow.h"
#include "Interactor.h"


MainWindow::MainWindow() {
    this -> setupUi();

    // A VTK renderer is created and bound with the qvtk widget
    this -> renderer = vtkSmartPointer<vtkRenderer>::New();
    this -> qvtkWidget -> GetRenderWindow() -> AddRenderer(renderer);
    this -> renderer -> SetGradientBackground (true);
    this -> renderer -> SetBackground (0.5, 0.5, 1);

}

void MainWindow::setupUi() {
    this -> resize(990, 583);

    // Central widget
    central_wrapping_widget = new QWidget();
    qvtkWidget = new QVTKWidget();
    layout_central = new QHBoxLayout();
    layout_central -> addWidget(qvtkWidget);
    central_wrapping_widget -> setLayout(layout_central);

    // Right Dock widget
    // dock_wrapping_widget = new QWidget();
    // menu_dock = new QDockWidget();
    // layout_dock_right = new QVBoxLayout();
    // dock_wrapping_widget -> setLayout(layout_dock_right);
    // menu_dock -> setWidget(dock_wrapping_widget);
    // addDockWidget(Qt::RightDockWidgetArea, menu_dock);

    //Top Menu Bar
    createActions();
    createMenus();

    this -> setCentralWidget(central_wrapping_widget);
    this -> setWindowTitle(QStringLiteral("PDART (WIP)"));


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
        vtkSmartPointer<vtkActorCollection> actor_coll = this -> renderer -> GetActors();
        int actor_count = actor_coll -> GetNumberOfItems();

        if (actor_count > 0) {
            for (int i = 0; i < actor_count; ++i)
                this -> renderer -> RemoveActor(actor_coll -> GetLastActor());
        }

        // The new actors are added to the scene
        this -> renderer -> AddActor(point_cloud_actor);
        this -> renderer -> AddActor(delaunay_actor);

        // The actions that were until now disabled are enabled
        this -> set_action_status(true, "&Shape Graphic Properties", "Shape Color");
        this -> set_action_status(true, "&Shape Graphic Properties", "Show Vertices");
        this -> set_action_status(true, "&Operation", "Select Points");
        this -> set_action_status(true, "&File", "&Save");


        // The interactor is set up and connected to the shape model being displayed
        vtkSmartPointer<InteractorStyle> style =
            vtkSmartPointer<InteractorStyle>::New();
        style -> SetPoints(read_polydata);
        style -> SetMainWindow(this);

        vtkSmartPointer<vtkAreaPicker> areaPicker =
            vtkSmartPointer<vtkAreaPicker>::New();
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();

        renderWindowInteractor -> SetPicker(areaPicker);
        renderWindowInteractor -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow());
        renderWindowInteractor -> SetInteractorStyle( style );



        this -> qvtkWidget -> GetRenderWindow() -> Render();



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
        shape_actor -> GetMapper()->ScalarVisibilityOff();
        shape_actor -> GetProperty()->SetColor(1, 1, 1);
        shape_actor -> GetProperty()->SetPointSize(3);


        // Remove any already existing actor from the rendering window
        vtkSmartPointer<vtkActorCollection> actor_coll = this -> renderer -> GetActors();
        int actor_count = actor_coll -> GetNumberOfItems();

        if (actor_count > 0) {
            for (int i = 0; i < actor_count; ++i)
                this -> renderer -> RemoveActor(actor_coll -> GetLastActor());
        }

        // The new actor is added to the scene
        this -> renderer -> AddActor(shape_actor);

        // The actions that were until now disabled are enabled
        this -> set_action_status(true, "&Shape Graphic Properties", "Shape Color");
        this -> set_action_status(true, "&Shape Graphic Properties", "Show Vertices");
        this -> set_action_status(true, "&Operation", "Select Points");
        this -> set_action_status(true, "&File", "&Save");


        // The interactor is set up and connected to the shape model being displayed
        vtkSmartPointer<InteractorStyle> style =
            vtkSmartPointer<InteractorStyle>::New();
        style -> SetPoints(read_polydata);
        style -> SetMainWindow(this);

        vtkSmartPointer<vtkAreaPicker> areaPicker =
            vtkSmartPointer<vtkAreaPicker>::New();
        vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
            vtkSmartPointer<vtkRenderWindowInteractor>::New();

        renderWindowInteractor -> SetPicker(areaPicker);
        renderWindowInteractor -> SetRenderWindow(this -> qvtkWidget -> GetRenderWindow());
        renderWindowInteractor -> SetInteractorStyle( style );

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

void MainWindow::setBackgroundColor() {
    QColorDialog * palette  = new QColorDialog();
    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> SetBackground (float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}

void MainWindow::setShapeColor() {
    QColorDialog * palette  = new QColorDialog();
    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> GetActors () -> GetLastActor()
        -> GetProperty() -> SetColor(float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}

void MainWindow::newShapeModel() {
    NewShapeModelDialog * newshapemodeldialog = new NewShapeModelDialog();

    if ( newshapemodeldialog -> exec() == QDialog::Accepted ) {
        int n_points = atoi(newshapemodeldialog -> input_box -> text().toStdString().c_str());
        vtkSmartPointer<vtkPointSource> new_pc = vtkSmartPointer<vtkPointSource>::New();
        new_pc -> SetNumberOfPoints(n_points);
        new_pc -> Update();
        vtkSmartPointer<vtkPolyData> new_polydata = new_pc -> GetOutput();
        this -> load_pc(new_polydata);
    }

}

void MainWindow::select() {
    selectorActive = true ;
    // TODO: bypass "r" key depression by simulating keyboard input
}

void MainWindow::save() {

    // The save path is queried
    QString fileName = QFileDialog::getSaveFileName(this,
                       tr("Save File"), "", tr("VTK PolyData file (*.vtp )"));
    if (!fileName.isEmpty()) {

        // the PolyData currently displayed is saved to a .vtp file
        // by means of the appropriate writer class
        vtkSmartPointer<vtkXMLPolyDataWriter> writer =
            vtkSmartPointer<vtkXMLPolyDataWriter>::New();
        writer -> SetFileName(fileName.toStdString().c_str());

        // Based on the insertion order in MainWindow::load(), the first actor in
        // the actor collection is the one we are interested in
        vtkSmartPointer<vtkActor> point_cloud_actor = vtkActor::SafeDownCast(this -> renderer -> GetActors() -> GetItemAsObject(0));
        vtkSmartPointer<vtkPolyDataMapper> pc_mapper = vtkPolyDataMapper::SafeDownCast(point_cloud_actor -> GetMapper());
        vtkSmartPointer<vtkDataSet> data_set = point_cloud_actor -> GetMapper() -> GetInput() ;

        // The polydata of interest is finally extracted
        vtkSmartPointer<vtkPolyData> polydata = vtkPolyData::SafeDownCast (data_set);


#if VTK_MAJOR_VERSION <= 5
        writer -> SetInput(polydata);
#else
        writer -> SetInputData(polydata);
#endif
        writer -> Write();

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
    connect(shapeColorAct, &QAction::triggered, this, &MainWindow::setShapeColor);
    shapeColorAct -> setDisabled(true);


    backgroundColorAct = new QAction(tr("Background Color"), this);
    backgroundColorAct -> setStatusTip(tr("Set the background color"));
    connect(backgroundColorAct, &QAction::triggered, this, &MainWindow::setBackgroundColor);

    newShapeModelAct = new QAction(tr("New"), this);
    newShapeModelAct -> setShortcuts(QKeySequence::New);
    newShapeModelAct -> setStatusTip(tr("Generate a new shape model"));
    connect(newShapeModelAct, &QAction::triggered, this, &MainWindow::newShapeModel);


    vertexVisibilityAct = new QAction(tr("Show Vertices"), this);
    vertexVisibilityAct -> setStatusTip(tr("Display or Hide Vertices"));
    vertexVisibilityAct -> setCheckable(true);
    vertexVisibilityAct -> setChecked(true);
    vertexVisibilityAct -> setDisabled(true);

    connect(vertexVisibilityAct, &QAction::triggered, this, &MainWindow::changeVertexVisibility);

}

void MainWindow::changeVertexVisibility() {
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
    fileMenu -> addAction(newShapeModelAct);
    fileMenu -> addSeparator();



    OperationMenu = menuBar()->addMenu(tr("&Operation"));
    OperationMenu -> addAction(selectPointAct);

    ViewMenu = menuBar()->addMenu(tr("&Shape Graphic Properties"));
    ViewMenu -> addAction(shapeColorAct);
    ViewMenu -> addAction(backgroundColorAct);
    ViewMenu -> addSeparator();
    ViewMenu -> addAction(vertexVisibilityAct);



}

vtkSmartPointer<vtkRenderer> MainWindow::getRenderer() {
    return this -> renderer;
}