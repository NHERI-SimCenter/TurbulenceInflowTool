#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "utilWindows/dialogabout.h"
#include "utilWindows/helpwindow.h"
#include <QGuiApplication>
#include <QScreen>
#include <QRect>

#include <QStandardItemModel>
#include <QDebug>

#include "customizeditemmodel.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->headerWidget->setHeadingText("Turbulent Inflow Adjustment Tool for OpenFOAM");

    standardModel = new CustomizedItemModel(); //QStandardItemModel ;
    QStandardItem *rootNode = standardModel->invisibleRootItem();

    //defining bunch of items for inclusion in model
    QStandardItem *sourceItem    = new QStandardItem("Source");
    QStandardItem *parameterItem = new QStandardItem("Parameters");
    QStandardItem *exportItem    = new QStandardItem("Export");

    //building up the hierarchy of the model
    rootNode->appendRow(sourceItem);
    rootNode->appendRow(parameterItem);
    rootNode->appendRow(exportItem);

    infoItemIdx = sourceItem->index();

    //register the model
    ui->treeView->setModel(standardModel);
    ui->treeView->expandAll();
    ui->treeView->setHeaderHidden(true);
    ui->treeView->setMinimumWidth(110);
    ui->treeView->setMaximumWidth(110);
    ui->treeView->setMinimumWidth(110);

    //Disable Edit for the TreeView
    ui->treeView->setEditTriggers(QTreeView::EditTrigger::NoEditTriggers);

    //
    // customize the apperance of the menu on the left
    //

    ui->treeView->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff ); // hide the horizontal scroll bar
    ui->treeView->setObjectName("treeViewOnTheLeft");
    ui->treeView->setIndentation(0);
    QFile file(":/styles/menuBar.qss");
    if (file.open(QFile::ReadOnly)) {
        ui->treeView->setStyleSheet(file.readAll());
        file.close();
    }
    else
    {
        qDebug() << "Open Style File Failed!";
    }

    //
    // set up so that a slection change triggers the selectionChanged slot
    //

    QItemSelectionModel *selectionModel= ui->treeView->selectionModel();
    connect(selectionModel,
            SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
            this,
            SLOT(selectionChangedSlot(const QItemSelection &, const QItemSelection &)));
    connect(ui->fileWidget, SIGNAL(hasValidSource(bool, QDir &)),
            ui->exportWidget, SLOT(setLocationAvailable(bool, QDir &)));
    connect(ui->exportWidget, SIGNAL(sendParameterMap()),
            ui->inflowWidget, SLOT(sendParameterMap()));
    connect(ui->inflowWidget, SIGNAL(parametersReady(QMap<QString, double> &)),
            ui->exportWidget, SLOT(setParameterMap(QMap<QString, double> &)));
    //
    // set active index
    //

    ui->treeView->setCurrentIndex(infoItemIdx);

    //
    // adjust size of application window to the available display
    //
    QRect rec = QGuiApplication::primaryScreen()->geometry();
    int height = this->height()<0.65*rec.height()?int(0.65*rec.height()):this->height();
    int width  = this->width()<0.65*rec.width()?int(0.65*rec.width()):this->width();
    this->resize(width, height);
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

void MainWindow::on_btn_selectSource_clicked()
{
    ui->inflowWidget->selectSourceLocation();
}

void
MainWindow::selectionChangedSlot(const QItemSelection & /*newSelection*/, const QItemSelection &/*oldSelection*/) {

    //get the text of the selected item
    const QModelIndex index = ui->treeView->selectionModel()->currentIndex();
    QString selectedText = index.data(Qt::DisplayRole).toString();

    if (selectedText == "Source")          { ui->theStackedWidget->setCurrentIndex(0); }
    else if (selectedText == "Parameters") { ui->theStackedWidget->setCurrentIndex(1); }
    else if (selectedText == "Export")     { ui->theStackedWidget->setCurrentIndex(2); }
    else {
        qWarning() << "Unknown page selected: " << selectedText;
        ui->theStackedWidget->setCurrentIndex(0);
    }
}

void MainWindow::on_action_Documentation_triggered()
{
    HelpWindow *help = new HelpWindow();
    help->show();
}
