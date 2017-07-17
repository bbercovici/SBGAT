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

    this -> setStatusBar(this -> status_bar);
    this -> statusBar() -> showMessage("Ready");

    this -> lateral_dockwidget -> setFeatures( QDockWidget::DockWidgetMovable );
    this -> lateral_dockwidget -> hide();

    createActions();
    createMenus();

    this -> setCentralWidget(qvtkWidget);
    this -> setWindowTitle(QStringLiteral("SBGAT (WIP)"));
    this -> addDockWidget(Qt::RightDockWidgetArea, this -> lateral_dockwidget);


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
            this -> create_vtkpolydata_from_shape_model(shape_model.get());



        }
    }




}

void Mainwindow::create_vtkpolydata_from_shape_model(SBGAT_CORE::ShapeModel * shape_model) {



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

    vtkSmartPointer<vtkActor> actor =
        vtkSmartPointer<vtkActor>::New();
    actor -> SetMapper(mapper);

    // Visualize
    this -> renderer -> AddActor(actor);

    this -> qvtkWidget -> GetRenderWindow() -> Render();







}

void Mainwindow::createMenus() {
    this -> fileMenu = menuBar()->addMenu(tr("&File"));
    this -> fileMenu -> addAction(this -> load_shape_model_action);

    this -> ViewMenu = menuBar() -> addMenu(tr("&Shape Graphic Properties"));
    this -> ViewMenu -> addAction(this -> set_background_color_action);
    this -> ViewMenu -> addSeparator();


}

vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}

