#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "utilWindows/dialogabout.h"
#include <QGuiApplication>
#include <QScreen>
#include <QRect>

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->header->setHeadingText("Turbulent Inflow Adjustment Tool for OpenFOAM");
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_action_Quit_triggered()
{
    this->close();
}

void MainWindow::on_action_About_triggered()
{
    DialogAbout *dlg = new DialogAbout("0.1");

    //
    // adjust size of application window to the available display
    //
    QRect rec = QGuiApplication::primaryScreen()->availableGeometry();
    int height = static_cast<int>(0.50*rec.height());
    int width  = static_cast<int>(0.50*rec.width());
    dlg->resize(width, height);

    dlg->exec();
    delete dlg;
}

void MainWindow::on_action_New_triggered()
{

}

void MainWindow::on_action_Open_triggered()
{

}

void MainWindow::on_action_Save_triggered()
{

}

void MainWindow::on_actionSave_As_triggered()
{

}
