#include "Mainwindow.hpp"

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
    this -> resize(990, 583);

    this -> lateral_dockwidget = new QDockWidget(this);
    this -> qvtkWidget = new QVTKWidget(this);
    this -> status_bar = new QStatusBar(this);
    this -> log_console = new QPlainTextEdit(this);
    this -> log_console -> setReadOnly(true);

    this -> setStatusBar(this -> status_bar);
    this -> statusBar() -> showMessage("Ready");

    // The lateral dockwidget is filled up
    this -> lateral_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );

    QWidget * lateral_dockwidget_container = new QWidget(this);
    QVBoxLayout * lateral_dockwidget_container_layout = new QVBoxLayout();

    lateral_dockwidget_container -> setLayout(lateral_dockwidget_container_layout);
    lateral_dockwidget_container_layout -> addWidget( this -> log_console );


    this -> lateral_dockwidget -> setWidget(lateral_dockwidget_container);
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> lateral_dockwidget);
    this -> lateral_dockwidget -> hide();

    // Central window
    this -> setCentralWidget(qvtkWidget);
    this -> setWindowTitle(QStringLiteral("SBGAT (WIP)"));


    // Actions and menus are created
    this -> createActions();
    this -> createMenus();



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

    this -> orientation_widget =
        vtkSmartPointer<vtkOrientationMarkerWidget>::New();
    orientation_widget -> SetOrientationMarker( axes );
    orientation_widget -> SetInteractor( this -> renderWindowInteractor );
    orientation_widget -> SetViewport( 0.0, 0.0, 0.2, 0.2 );
    orientation_widget -> SetEnabled( 1 );
    orientation_widget -> InteractiveOff();


    this -> show();
    this -> qvtkWidget -> GetRenderWindow() -> Render();

}

vtkSmartPointer<vtkRenderWindowInteractor> Mainwindow::get_render_window_interactor() {
    return this -> renderWindowInteractor;
}

void Mainwindow::set_action_status(bool enabled, QAction * action) {
    action -> setEnabled(enabled);
}


void Mainwindow::set_background_color() {
    QColorDialog *  palette = new QColorDialog(this);

    QColor qcolor =  palette -> getColor();
    if (qcolor.isValid()) {
        this -> renderer -> SetBackground (float(qcolor.red()) / 255, float(qcolor.green()) / 255, float(qcolor.blue()) / 255);
        this -> qvtkWidget -> GetRenderWindow() -> Render();
    }

}


void Mainwindow::createActions() {

    this -> set_background_color_action = new QAction(tr("Background Color"), this);
    this -> set_background_color_action -> setStatusTip(tr("Set the background color"));
    connect(this -> set_background_color_action, &QAction::triggered, this, &Mainwindow::set_background_color);


    this -> load_shape_model_action = new QAction(tr("Load shape model"), this);
    this -> load_shape_model_action -> setStatusTip(tr("Load obj file holding the facet/vertex description of a shape of interest"));
    connect(this -> load_shape_model_action, &QAction::triggered, this, &Mainwindow::load_shape_model);


    this -> show_lateral_dockwidget_action = new QAction(tr("Show lateral widget"), this);
    this -> show_lateral_dockwidget_action -> setStatusTip(tr("Shows/hides lateral widget holding shape model information"));
    connect(this -> show_lateral_dockwidget_action, &QAction::triggered, this, &Mainwindow::show_lateral_dockwidget);


    this -> clear_console_action = new QAction(tr("Clear log console"), this);
    this -> clear_console_action -> setStatusTip(tr("Clears the log console"));
    connect(this -> clear_console_action, &QAction::triggered, this, &Mainwindow::clear_console);


    this -> save_console_action = new QAction(tr("Save log console to file"), this);
    this -> save_console_action -> setStatusTip(tr("Saves log console to a file"));
    connect(this -> save_console_action, &QAction::triggered, this, &Mainwindow::save_console);


    this -> print_inertia_action = new QAction(tr("Get inertia tensor"), this);
    this -> print_inertia_action -> setStatusTip(tr("Print inertia tensor to the log console"));
    connect(this -> print_inertia_action, &QAction::triggered, this, &Mainwindow::print_inertia);

    this -> compute_pgm_acceleration_action = new QAction(tr("Compute PGM acceleration"), this);
    this -> compute_pgm_acceleration_action -> setStatusTip(tr("Compute PGM acceleration at a point whose coordinates are expressed in the shape's body frame"));
    connect(this -> compute_pgm_acceleration_action, &QAction::triggered, this, &Mainwindow::compute_pgm_acceleration);



    this -> compute_global_pgm_acceleration_action = new QAction(tr("Compute global PGM accelerations"), this);
    this -> compute_global_pgm_acceleration_action -> setStatusTip(tr("Compute PGM accelerations over the entire shape model"));
    connect(this -> compute_global_pgm_acceleration_action, &QAction::triggered, this, &Mainwindow::compute_global_pgm_acceleration);


    this -> compute_gravity_slopes_action = new QAction(tr("Compute gravity slopes"), this);
    this -> compute_gravity_slopes_action -> setStatusTip(tr("Compute PGM accelerations over the entire shape model"));
    connect(this -> compute_gravity_slopes_action, &QAction::triggered, this, &Mainwindow::compute_gravity_slopes);


    this -> toggle_grav_slopes_visibility_action = new QAction(tr("Show gravity slopes"), this);
    this -> toggle_grav_slopes_visibility_action -> setStatusTip(tr("Display gravity slopes along with colorbar"));
    connect(this -> toggle_grav_slopes_visibility_action, &QAction::triggered, this, &Mainwindow::toggle_grav_slopes_visibility);

}


void Mainwindow::toggle_grav_slopes_visibility() {

    vtkSmartPointer<vtkPolyData> active_polydata = this -> polydatas.begin() -> second;
    vtkSmartPointer<vtkPolyDataMapper> active_mapper = this -> mappers.begin() -> second;



    if (active_mapper -> GetScalarVisibility()) {
        active_mapper -> ScalarVisibilityOff();
        this -> toggle_grav_slopes_visibility_action -> setText("Show gravity slopes");
        if (this -> renderer -> GetActors2D() -> GetNumberOfItems () > 0) {
            this -> renderer -> RemoveActor2D(this -> renderer -> GetActors2D() -> GetLastActor2D());
        }

    }
    else {
        active_mapper -> ScalarVisibilityOn();
        active_mapper -> SetScalarModeToUseCellData();



        double range[2] ;
        active_polydata -> GetCellData() -> GetScalars() -> GetRange(range);

        active_mapper -> SetColorModeToMapScalars();

        active_mapper -> SetScalarRange(range[0], range[1]);

        vtkSmartPointer<vtkScalarBarActor> scalarBar =
            vtkSmartPointer<vtkScalarBarActor>::New();
        scalarBar -> SetLookupTable(active_mapper -> GetLookupTable());
        scalarBar -> SetTitle("Gravity slopes (deg)");
        scalarBar -> SetNumberOfLabels(4);



        this -> renderer -> AddActor2D(scalarBar);

        this -> toggle_grav_slopes_visibility_action -> setText("Hide gravity slopes");


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


void Mainwindow::show_lateral_dockwidget() {

    if (this -> lateral_dockwidget -> isVisible()) {
        this -> lateral_dockwidget -> hide();
        this -> show_lateral_dockwidget_action -> setText(QString::fromStdString("Show lateral widget"));
    }
    else {
        this -> lateral_dockwidget -> show();
        this -> show_lateral_dockwidget_action -> setText(QString::fromStdString("Hide lateral widget"));

    }


}


void Mainwindow::load_shape_model() {

    QString fileName = QFileDialog::getOpenFileName(this,
                       tr("Open Shape Model"), "~/", tr("Wavefront file (*.obj)"));

    if (fileName.isEmpty() == false) {

        bool ok;
        double scaling_factor = QInputDialog::getDouble(this,
                                "Scaling factor", "Enter scaling factor :", 1, 1e-6, 1e6,
                                5, &ok);
        if (ok) {

            std::chrono::time_point<std::chrono::system_clock> start, end;

            start = std::chrono::system_clock::now();
            SBGAT_CORE::ShapeModelImporter shape_io(fileName.toStdString(),
                                                    scaling_factor, true);

            std::shared_ptr<SBGAT_CORE::ShapeModel>  shape_model = std::make_shared<SBGAT_CORE::ShapeModel>("BF", this -> frame_graph.get());

            shape_io.load_shape_model(shape_model.get());

            // The name of the shape model is extracted from the path
            int dot_index = fileName.lastIndexOf(".");
            int slash_index = fileName.lastIndexOf("/");
            std::string name = (fileName.toStdString()).substr(slash_index + 1 , dot_index - slash_index - 1);

            // The shape loaded is stored
            this -> shape_models[name] = shape_model;

            // A VTK Polydata is created from the loaded shape model
            // and displayed on the QVTKWidget
            this -> create_vtkpolydata_from_shape_model(shape_model.get(), name);
            end = std::chrono::system_clock::now();

            std::chrono::duration<double> elapsed_seconds = end - start;

            // The status bar is updated to show the significant figures of the current shape model
            std::string message("Facets : " + std::to_string(shape_model -> get_NFacets()) + " Vertices: " + std::to_string(shape_model -> get_NVertices())
                                + " Edges: " + std::to_string(shape_model -> get_NEdges()));

            this -> statusBar() -> showMessage(QString::fromStdString(message));
            if (this -> lateral_dockwidget -> isVisible() == false) {
                this -> show_lateral_dockwidget();

            }

            // The log console displays the name and content of the loaded shape model
            this -> log_console -> appendPlainText(QString::fromStdString("- Loading shape model from ") + fileName);
            this -> log_console -> appendPlainText(QString::fromStdString("- Loading completed in ")
                                                   + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));

        }
    }

}

void Mainwindow::print_inertia() const {

    this -> log_console -> appendPlainText(QString::fromStdString("- Dimensionless principal inertia matrix: "));
    std::stringstream ss;


    this -> shape_models.begin() -> second -> get_inertia().print(ss);

    this -> log_console -> appendPlainText(QString::fromStdString(ss.str()));


}






void Mainwindow::create_vtkpolydata_from_shape_model(SBGAT_CORE::ShapeModel * shape_model, std::string name) {


    // Add the polygon to a list of polygons
    vtkSmartPointer<vtkCellArray> polygons =
        vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkPoints> points =
        vtkSmartPointer<vtkPoints>::New();



    for (unsigned int facet_index = 0; facet_index < shape_model -> get_NFacets(); ++facet_index) {

        auto vertices = shape_model -> get_facets() -> at(facet_index) -> get_vertices() ;

        double * p0 = vertices -> at(0) -> get_coordinates() -> colptr(0);
        double * p1 = vertices -> at(1) -> get_coordinates() -> colptr(0);
        double * p2 = vertices -> at(2) -> get_coordinates() -> colptr(0);

        points -> InsertNextPoint(p0);
        points -> InsertNextPoint(p1);
        points -> InsertNextPoint(p2);

        // Create the polygon
        vtkSmartPointer<vtkPolygon> polygon =
            vtkSmartPointer<vtkPolygon>::New();
        polygon -> GetPointIds() -> SetNumberOfIds(3); //make a triangle
        polygon -> GetPointIds() -> SetId(0, 3 * facet_index);
        polygon -> GetPointIds() -> SetId(1, 3 * facet_index + 1);
        polygon -> GetPointIds() -> SetId(2, 3 * facet_index + 2);

        polygons -> InsertNextCell(polygon);


    }




    // Create a PolyData
    vtkSmartPointer<vtkPolyData> polygonPolyData =
        vtkSmartPointer<vtkPolyData>::New();
    polygonPolyData -> SetPoints(points);
    polygonPolyData -> SetPolys(polygons);




    // Create a mapper and actor
    vtkSmartPointer<vtkPolyDataMapper> mapper =
        vtkSmartPointer<vtkPolyDataMapper>::New();

    mapper -> SetInputData(polygonPolyData);
    mapper -> ScalarVisibilityOff();

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor -> SetMapper(mapper);

    // Visualize
    this -> renderer -> AddActor(actor);


    this -> qvtkWidget -> GetRenderWindow() -> Render();

    // Store
    this -> polydatas[name] = polygonPolyData;
    this -> mappers[name] = mapper;
    this -> actors[name] = actor;




}








void Mainwindow::compute_global_pgm_acceleration() {

    SBGAT_CORE::DynamicAnalyses dynas(this -> shape_models.begin() -> second.get());

    bool ok_density = false;


    double density = QInputDialog::getDouble(this,
                     "Global Polyhedron Gravity Model Acceleration", "Density (kg/m^3) :", 2000, 0, 1e9,
                     5, &ok_density);

    if (ok_density) {



        this -> log_console -> appendPlainText(QString::fromStdString("- Computing global accelerations..."));


        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();
        dynas.compute_pgm_accelerations(density);
        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;

        this -> log_console -> appendPlainText(QString::fromStdString("- Done computing in ")
                                               + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));



    }

}

void Mainwindow::compute_gravity_slopes() {

    SBGAT_CORE::DynamicAnalyses dynas(this -> shape_models.begin() -> second.get());

    bool ok_spin_axis = true;
    bool correct_format = false;
    double spin_rate;
    arma::vec angles_arma(3);
    bool ok_spin_rate = false;
    arma::vec spin_axis = {0, 0, 1};


    while (ok_spin_axis == true && correct_format == false) {

        QString coords = QInputDialog::getText(this,
                                               tr("Gravity slopes"),
                                               tr("313 Euler angles (deg) directing spin axis ((0,0,0) : along z) :"),
                                               QLineEdit::Normal,
                                               QString::fromStdString("Omega,i,omega") ,
                                               &ok_spin_axis);

        QStringList angles_split = coords.split(",");

        if (angles_split.count() != 3) {
            correct_format = false;
            continue;
        }

        // Matching doubles
        QRegExp re("[-+]?[0-9]*\.?[0-9]+");

        if (re.exactMatch(angles_split.at(0)) && re.exactMatch(angles_split.at(1)) && re.exactMatch(angles_split.at(2))) {
            angles_arma(0) = angles_split.at(0).toDouble();
            angles_arma(1) = angles_split.at(1).toDouble();
            angles_arma(2) = angles_split.at(2).toDouble();

            correct_format = true;


        }
        else {
            correct_format = false;
            continue;
        }

    }



    if (ok_spin_axis) {
        spin_axis = RBK::euler313_to_dcm(angles_arma).t() * spin_axis;

        spin_rate = QInputDialog::getDouble(this,
                                            "Gravity slopes", "Angular velocity (rad/s) :", 0, -1e9, 1e9,
                                            5, &ok_spin_rate);

        if (ok_spin_rate ) {

            this -> log_console -> appendPlainText(QString::fromStdString("- Computing gravity slopes..."));


            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();
            dynas.compute_gravity_slopes(spin_axis, spin_rate);

            this -> update_vtk_slopes();

            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;





            this -> log_console -> appendPlainText(QString::fromStdString("- Done computing in ")
                                                   + QString::number(elapsed_seconds.count()) +  QString::fromStdString(" s"));


        }

    }
}


void Mainwindow::update_vtk_slopes() {

    vtkSmartPointer<vtkPolyData> active_polydata = this -> polydatas.begin() -> second;
    vtkSmartPointer<vtkPolyDataMapper> active_mapper = this -> mappers.begin() -> second;
    std::shared_ptr<SBGAT_CORE::ShapeModel> active_shape_model = this -> shape_models.begin() -> second;

    vtkSmartPointer<vtkDoubleArray> slope_data =
        vtkSmartPointer<vtkDoubleArray>::New();
    slope_data -> SetNumberOfValues(active_shape_model -> get_NFacets());
    slope_data -> SetName("SlopeData");

    for (unsigned int facet_index = 0; facet_index < active_shape_model -> get_NFacets(); ++facet_index) {

        SBGAT_CORE::Facet * facet =  active_shape_model -> get_facets() -> at(facet_index);

        slope_data -> SetValue(facet_index, facet -> get_facet_results() -> get_grav_slope());
    }

    // The array is added to the polydata
    active_polydata -> GetCellData() -> SetScalars(slope_data);
    active_polydata -> Modified();

}


void Mainwindow::compute_pgm_acceleration() {


    SBGAT_CORE::DynamicAnalyses dynas(this -> shape_models.begin() -> second.get());

    bool ok_coords = true;
    bool correct_format = false;
    double point[3];
    bool ok_density = false;
    double density;

    while (ok_coords == true && correct_format == false) {

        QString coords = QInputDialog::getText(this,
                                               tr("Polyhedron Gravity Model Acceleration"),
                                               tr("Body-fixed frames coordinates (m) :"),
                                               QLineEdit::Normal,
                                               QString::fromStdString("x,y,z") ,
                                               &ok_coords);

        QStringList coords_split = coords.split(",");

        if (coords_split.count() != 3) {
            correct_format = false;
            continue;
        }

        // Matching doubles
        QRegExp re("[-+]?[0-9]*\.?[0-9]+");

        if (re.exactMatch(coords_split.at(0)) && re.exactMatch(coords_split.at(1)) && re.exactMatch(coords_split.at(2))) {
            point[0] = coords_split.at(0).toDouble();
            point[1] = coords_split.at(1).toDouble();
            point[2] = coords_split.at(2).toDouble();

            correct_format = true;


        }
        else {
            correct_format = false;
            continue;
        }

    }



    if (ok_coords) {
        density = QInputDialog::getDouble(this,
                                          "Polyhedron Gravity Model Acceleration", "Density (kg/m^3) :", 2000, 0, 1e9,
                                          5, &ok_density);


        if (ok_density) {

            // The PGM acceleration is computed at the provided point
            arma::vec coords_arma = {point[0], point[1], point[2]};
            std::stringstream ss_coords;

            arma::vec acc = dynas.pgm_acceleration(point , density);
            std::stringstream ss_acc;

            acc.print(ss_acc);
            coords_arma.print(ss_coords);

            this -> log_console -> appendPlainText(QString::fromStdString("\n- At body-fixed coordinates (m) : "));
            this -> log_console -> appendPlainText(QString::fromStdString(ss_coords.str()));
            this -> log_console -> appendPlainText(QString::fromStdString("- PGM acceleration (m/s^2): "));
            this -> log_console -> appendPlainText(QString::fromStdString(ss_acc.str()));

        }
    }




}


void Mainwindow::createMenus() {
    this -> FileMenu = menuBar()->addMenu(tr("&File"));
    this -> FileMenu -> addAction(this -> load_shape_model_action);

    this -> ViewMenu = menuBar() -> addMenu(tr("&View"));
    this -> ViewMenu -> addAction(this -> set_background_color_action);
    this -> ViewMenu -> addAction(this -> show_lateral_dockwidget_action);
    this -> ViewMenu -> addSeparator();

    this -> ShapeMenu = menuBar() -> addMenu(tr("&Shape Model"));
    this -> ShapeMenu -> addAction(this -> print_inertia_action);

    this -> DynamicAnalysesMenu = menuBar() -> addMenu(tr("&Dynamic Analyses"));
    this -> DynamicAnalysesMenu -> addAction(this -> compute_pgm_acceleration_action);
    this -> DynamicAnalysesMenu -> addSeparator();
    this -> DynamicAnalysesMenu -> addAction(this -> compute_global_pgm_acceleration_action);
    this -> DynamicAnalysesMenu -> addAction(this -> compute_gravity_slopes_action);


    this -> ResultsMenu = menuBar() -> addMenu(tr("&Results"));
    this -> ResultsMenu -> addAction(this -> toggle_grav_slopes_visibility_action);

    this -> ConsoleMenu = menuBar() -> addMenu(tr("&Log Console"));
    this -> ConsoleMenu -> addAction(this -> clear_console_action);
    this -> ConsoleMenu -> addAction(this -> save_console_action);







}

vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}

