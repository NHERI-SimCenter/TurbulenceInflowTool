#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "utilWindows/dialogabout.h"
#include "utilWindows/helpwindow.h"
#include <QGuiApplication>
#include <QScreen>
#include <QRect>
#include <QMessageBox>
#include <QDesktopServices>
#include <QUrl>
#include <QFileDialog>
#include <QJsonObject>
#include <QJsonDocument>
#include <QDateTime>
#include <QTextBrowser>

#include <QStandardItemModel>
#include <QDebug>

#include "customizeditemmodel.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->headerWidget->setHeadingText("Turbulence Inflow Tool");

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

#ifdef Q_OS_WIN
#define COLWIDTH 140
#else
#define COLWIDTH 110
#endif

    //register the model
    ui->treeView->setModel(standardModel);
    ui->treeView->expandAll();
    //ui->treeView->setHeaderHidden(true);
    //ui->treeView->header()->setStretchLastSection(true);
    ui->treeView->setMinimumWidth(COLWIDTH);
    ui->treeView->setMaximumWidth(COLWIDTH);
    ui->treeView->setColumnWidth(0,COLWIDTH);
    ui->treeView->setIconSize(QSize(0,0));

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

    connect(selectionModel,   SIGNAL(selectionChanged(const QItemSelection &, const QItemSelection &)),
            this,             SLOT(selectionChangedSlot(const QItemSelection &, const QItemSelection &)));
    connect(ui->fileWidget,   SIGNAL(hasValidSource(bool, QDir &)),
            ui->exportWidget, SLOT(setLocationAvailable(bool, QDir &)));
    connect(ui->exportWidget, SIGNAL(sendParameterMap()),
            ui->inflowWidget, SLOT(sendParameterMap()));
    connect(ui->inflowWidget, SIGNAL(parametersReady(QMap<QString, double> &)),
            ui->exportWidget, SLOT(setParameterMap(QMap<QString, double> &)));
    connect(ui->fileWidget,   SIGNAL(hasValidSource(bool, QDir &)),
            this,             SLOT(fetchUFileData(bool, QDir &)));
    connect(ui->fileWidget,   SIGNAL(boundarySelection(int)),
            ui->exportWidget, SLOT(setBoundarySelection(int)));
    connect(ui->exportWidget, SIGNAL(boundarySelection(int)),
            ui->fileWidget,   SLOT(setBoundarySelection(int)));
    connect(ui->fileWidget,   SIGNAL(sendModel(QStandardItemModel *)),
            ui->exportWidget, SLOT(setModel(QStandardItemModel *)));

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

    // setting text
    versionText = "Turbulence Inflow Tool - Version " + QString(APP_VERSION);
    citeText = "Jiawei Wan, Peter Mackenzie-Helnwein, and Frank McKenna. (2019, September 26). NHERI-SimCenter/TurbulentInflowTool: Versions 1.0.2 (Version v1.0.2). Zenodo. http://doi.org/10.5281/zenodo.3462805";

    manualURL = "https://www.designsafe-ci.org/data/browser/public/designsafe.storage.community//SimCenter/Software/TurbulentInflowTool/";
    feedbackURL = "https://docs.google.com/forms/d/e/1FAIpQLSfh20kBxDmvmHgz9uFwhkospGLCeazZzL770A2GuYZ2KgBZBA/viewform";
    featureRequestURL = "https://docs.google.com/forms/d/e/1FAIpQLScTLkSwDjPNzH8wx8KxkyhoIT7AI9KZ16Wg9TuW1GOhSYFOag/viewform";
    copyrightText = QString("\
                            <p>\
                            The source code is licensed under a BSD 2-Clause License:<p>\
                            \"Copyright (c) 2017-2018, The Regents of the University of California (Regents).\"\
                            All rights reserved.<p>\
                            <p>\
                            Redistribution and use in source and binary forms, with or without \
                            modification, are permitted provided that the following conditions are met:\
                            <p>\
                            1. Redistributions of source code must retain the above copyright notice, this\
                            list of conditions and the following disclaimer.\
                            \
                            \
                            2. Redistributions in binary form must reproduce the above copyright notice,\
                            this list of conditions and the following disclaimer in the documentation\
                            and/or other materials provided with the distribution.\
                            <p>\
                            THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS \"AS IS\" AND\
                            ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED\
                            WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE\
                            DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR\
                            ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES\
                            (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;\
            LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND\
            ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT\
            (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS\
            SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\
            <p>\
            The views and conclusions contained in the software and documentation are those\
            of the authors and should not be interpreted as representing official policies,\
            either expressed or implied, of the FreeBSD Project.\
            <p>\
            REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, \
            THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.\
            THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS \
            PROVIDED \"AS IS\". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT,\
            UPDATES, ENHANCEMENTS, OR MODIFICATIONS.\
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            The compiled binary form of this application is licensed under a GPL Version 3 license.\
            The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
            included in the packaging of this application. \
            <p>\
            ------------------------------------------------------------------------------------\
            <p>\
            This software makes use of the QT packages (unmodified): core, gui, widgets and network\
                                                                     <p>\
                                                                     QT is copyright \"The Qt Company Ltd&quot; and licensed under the GNU Lesser General \
                                                                     Public License (version 3) which references the GNU General Public License (version 3)\
            <p>\
            The licenses are as published by the Free Software Foundation and appearing in the LICENSE file\
            included in the packaging of this application. \
            <p>\
      ");

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
    DialogAbout *dlg = new DialogAbout(versionText);

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
    ui->inflowWidget->reset();
}

void MainWindow::on_action_Open_triggered()
{
    /* identify filename and location for loading */

    QString theFolder = QDir::homePath();
    QString theFilter = "Json file (*.json)";
    QFileDialog dlg;

    QString filename = dlg.getOpenFileName(this, "Load file", theFolder, theFilter);

    qWarning() << filename;

    /* load JSON object from file */
    QFile loadFile;
    loadFile.setFileName(filename);

    if (!loadFile.open(QIODevice::ReadOnly | QIODevice::Text)) {
        QMessageBox msg(QMessageBox::Information, "Info", "Could not open file.");
        msg.exec();
        return;
    }

    QString theFile = loadFile.readAll();
    loadFile.close();

    bool sourceValid = false;

    QJsonDocument infoDoc = QJsonDocument::fromJson(theFile.toUtf8());

    /* start a JSON object to represent the system */
    QJsonObject json = infoDoc.object();

    QString creator;
    creator  = json["creator"].toString();

    if (creator == "TurbulentInflowTool")  sourceValid = true;  // v1.0.0 and 1.0.1
    if (creator == "TurbulenceInflowTool") sourceValid = true;  // >= v1.0.2

    QString version;
    version  = json["version"].toString();

    bool versionValid = false;
    if (version.startsWith("1.") ) versionValid = true;
    //if (version.startsWith("2.0") ) versionValid = true;
    //if (version.startsWith("2.1") ) versionValid = true;

    if (sourceValid && versionValid) {

        QString username;
        username = json["username"].toString();
        QString author;
        author   = json["author"].toString();
        QString filedate;
        filedate = json["date"].toString();

        /* read parameter information and update UI */

        if (version.startsWith("1.0"))
        {
            if (json.contains("parameters")) {
                QJsonObject params = json["parameters"].toObject();
                ui->inflowWidget->inputFromJSON(params);
            }
        }
        else
        {
            QMessageBox msg(QMessageBox::Information, "Info", "Could not read file: invalid Version");
            msg.exec();
        }
    }
    else
    {
        QMessageBox msg(QMessageBox::Information, "Info", "Not a valid model file.");
        msg.exec();
    }

}

void MainWindow::on_action_Save_triggered()
{
    /* identify filename and location for saving */

    QString path = QDir::homePath();
    QDir d;
    d.setPath(path);
    QString filename = d.filePath("TurbulenceInflowTool.json");
    QString theFilter = "Json file (*.json)";
    QFileDialog dlg;

    filename = dlg.getSaveFileName(this, "Save file", filename, theFilter );

    // check if cancelled
    if (!filename.isEmpty())
    {
        /* start a JSON object to represent the system */
        QJsonObject *json = new QJsonObject();

        json->insert("creator", QString("TurbulenceInflowTool"));
        json->insert("version", QString(APP_VERSION));
#ifdef Q_OS_WIN
        QString username = qgetenv("USERNAME");
#else
        QString username = qgetenv("USER");
#endif
        json->insert("author", username);
        json->insert("date", QDateTime::currentDateTime().toString());

        QJsonObject params;

        ui->inflowWidget->outputToJSON(params);

        json->insert("parameters", params);

        QJsonDocument infoDoc = QJsonDocument(*json);

        /* write JSON object to file */

        QFile saveFile( filename );

        if (saveFile.open(QIODevice::WriteOnly)) {

            saveFile.write( infoDoc.toJson() );
            saveFile.close();
        }
        else
        {
            QMessageBox msg(QMessageBox::Information, "Info", "Could not save to file.");
            msg.exec();
        }

        // clean up
        delete json;
    }
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
    QDesktopServices::openUrl(QUrl(manualURL, QUrl::TolerantMode));
    //QDesktopServices::openUrl(QUrl("https://www.designsafe-ci.org/help/new-ticket/", QUrl::TolerantMode));

    //HelpWindow *help = new HelpWindow();
    //help->show();
}


void MainWindow::fetchUFileData(bool, QDir &)
{
    if (ui->fileWidget->fetchUFileData(head, tail, data))
    {
        ui->exportWidget->setUFileData(head, tail, data);
    }
}

void MainWindow::on_actionLicense_triggered()
{
    QDialog msgBox;
    QTextBrowser *theBrowser = new QTextBrowser();
    theBrowser->setSizePolicy(QSizePolicy::Expanding, QSizePolicy::Expanding);
    theBrowser->setText(copyrightText);
    QGridLayout *layout = new QGridLayout();
    layout->addWidget(theBrowser);
    msgBox.setLayout(layout);

    QRect rec = QGuiApplication::primaryScreen()->geometry();
    int height = this->height()<0.80*rec.height()?int(0.80*rec.height()):this->height();
    int width  = this->width()<0.65*rec.width()?int(0.65*rec.width()):this->width();
    msgBox.resize(width, height);

    msgBox.exec();
}

void MainWindow::on_actionHow_to_cite_triggered()
{
    QMessageBox msgBox;
    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(citeText);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());

    msgBox.exec();
}

void MainWindow::on_actionProvide_feeback_triggered()
{
    QDesktopServices::openUrl(QUrl(feedbackURL, QUrl::TolerantMode));
    //QDesktopServices::openUrl(QUrl("https://www.designsafe-ci.org/help/new-ticket/", QUrl::TolerantMode));
}

void MainWindow::on_actionSubmit_Feature_Request_triggered()
{
    QDesktopServices::openUrl(QUrl(featureRequestURL, QUrl::TolerantMode));
    //QDesktopServices::openUrl(QUrl("https://www.designsafe-ci.org/help/new-ticket/", QUrl::TolerantMode));
}

void MainWindow::on_action_Version_triggered()
{
    QMessageBox msgBox;
    QSpacerItem *theSpacer = new QSpacerItem(700, 0, QSizePolicy::Minimum, QSizePolicy::Expanding);
    msgBox.setText(versionText);
    QGridLayout *layout = (QGridLayout*)msgBox.layout();
    layout->addItem(theSpacer, layout->rowCount(),0,1,layout->columnCount());
    msgBox.exec();
}

