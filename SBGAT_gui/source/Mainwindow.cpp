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


}


void Mainwindow::createMenus() {
    this -> fileMenu = menuBar()->addMenu(tr("&File"));

    this -> ViewMenu = menuBar() -> addMenu(tr("&Shape Graphic Properties"));
    this -> ViewMenu -> addAction(this -> set_background_color_action);
    this -> ViewMenu -> addSeparator();


}

vtkSmartPointer<vtkRenderer> Mainwindow::get_renderer() {
    return this -> renderer;
}

