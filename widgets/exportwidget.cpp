#include "exportwidget.h"
#include "ui_exportwidget.h"

#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QMap>
#include <QDir>

#include <QDebug>

ExportWidget::ExportWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::ExportWidget)
{
    ui->setupUi(this);
    theParameters.clear();
    hasParameters = false;
}

ExportWidget::~ExportWidget()
{
    delete ui;
}

void ExportWidget::setLocationAvailable(bool status, QDir &loc)
{
    qDebug() << "exportWidget: signal received";

    if (status) {
        hasLocation = true;
        oldLocation = loc;
        newLocation = loc;
    }
    else {
        hasLocation = false;
        oldLocation = QDir::homePath();
        newLocation = QDir::homePath();
    }
}


void ExportWidget::setParameterMap(QMap<QString, double> &map)
{
    theParameters = map;
    hasParameters = true;
}

void ExportWidget::exportInflowParameterFile(QString fileName)
{
    hasParameters = false;

    // requests parameters to be sent
    emit sendParameterMap();

    // wait for parameters to arrive
    int i = 0;
    while (!hasParameters) { i++; }

    qDebug() << "Had to wait for " << i << "cycles";

    QString profile;

    switch (int(theParameters.value("profile")))
    {
    case 0: { profile="linear"; break; }
    case 1: { profile="exponential"; break; }
    case 2: { profile="linear"; break; }
    case 3: { profile="exponential"; break; }
    default: { profile="linear"; break; }
    }

    QFile theFile(fileName);
    if (theFile.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream out(&theFile);
        out << "Result: " << qSetFieldWidth(10) << left << 3.14 << 2.7;
        // writes "Result: 3.14      2.7       "

        out << "/*--------------------------------*- C++ -*----------------------------------*\\";
        out << "  =========                 |";
        out << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox";
        out << "   \\    /   O peration     | Website:  https://openfoam.org";
        out << "    \\  /    A nd           | Version:  6";
        out << "     \\/     M anipulation  |";
        out << "\\*---------------------------------------------------------------------------*/";
        out << "FoamFile";
        out << "{";
        out << "    version     2.0;";
        out << "    format      ascii;";
        out << "    class       dictionary;";
        out << "    location    \"constant\";";
        out << "    object      inflowProperties;";
        out << "}";
        out << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //";
        out << "";
        out << "        // mean velocity";
        out << "        U_Profile";
        out << "        {";

        out << "                        referenceValue    ("
            << theParameters.value("vel0")
            << theParameters.value("vel1")
            << theParameters.value("vel2")
            << ");";

        out << "                        profile           " << profile << ";";

        out << "                        direction         ("
            << theParameters.value("dir0")
            << theParameters.value("dir1")
            << theParameters.value("dir2")
            << ");";

        out << "                        origin            ("
            << theParameters.value("x0")
            << theParameters.value("x1")
            << theParameters.value("x2")
            << ");";

        out << "                        referencePoint    ("
            << theParameters.value("referencePoint0")
            << theParameters.value("referencePoint1")
            << theParameters.value("referencePoint2")
            << ");";

        out << "                        alpha             ("
            << theParameters.value("alpha0")
            << theParameters.value("alpha1")
            << theParameters.value("alpha2")
            << ");";

        out << "        }";
        out << "";
        out << "";
        out << "        // turbulent intensity (symmTensorField)";
        out << "        I_Profile";
        out << "        {";

        out << "                        profile           " << profile << ";";

        out << "                        referenceValue    ("
            << theParameters.value("phi00")
            << theParameters.value("phi10")
            << theParameters.value("phi20")
            << theParameters.value("phi11")
            << theParameters.value("phi21")
            << theParameters.value("phi22")
            << ");";

        out << "                        direction         ("
            << theParameters.value("dir0")
            << theParameters.value("dir1")
            << theParameters.value("dir2")
            << ");";

        out << "                        origin            ("
            << theParameters.value("x0")
            << theParameters.value("x1")
            << theParameters.value("x2")
            << ");";

        out << "                        referencePoint    ("
            << theParameters.value("referencePoint0")
            << theParameters.value("referencePoint1")
            << theParameters.value("referencePoint2")
            << ");";

        out << "                        alpha             ("
            << theParameters.value("alpha0")
            << theParameters.value("alpha1")
            << theParameters.value("alpha2")
            << ");";

        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for u component";
        out << "        Lu_Profile";
        out << "        {";

        out << "                        profile           " << profile << ";";

        out << "                        referenceValue    ("
            << theParameters.value("Lu0")
            << theParameters.value("Lu1")
            << theParameters.value("Lu2")
            << ");";

        out << "                        direction         ("
            << theParameters.value("dir0")
            << theParameters.value("dir1")
            << theParameters.value("dir2")
            << ");";

        out << "                        origin            ("
            << theParameters.value("x0")
            << theParameters.value("x1")
            << theParameters.value("x2")
            << ");";

        out << "                        referencePoint    ("
            << theParameters.value("referencePoint0")
            << theParameters.value("referencePoint1")
            << theParameters.value("referencePoint2")
            << ");";

        out << "                        alpha             ("
            << theParameters.value("alpha0")
            << theParameters.value("alpha1")
            << theParameters.value("alpha2")
            << ");";

        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for v component";
        out << "        Lv_Profile";
        out << "        {";

        out << "                        profile           " << profile << ";";

        out << "                        referenceValue    ("
            << theParameters.value("Lv0")
            << theParameters.value("Lv1")
            << theParameters.value("Lv2")
            << ");";

        out << "                        direction         ("
            << theParameters.value("dir0")
            << theParameters.value("dir1")
            << theParameters.value("dir2")
            << ");";

        out << "                        origin            ("
            << theParameters.value("x0")
            << theParameters.value("x1")
            << theParameters.value("x2")
            << ");";

        out << "                        referencePoint    ("
            << theParameters.value("referencePoint0")
            << theParameters.value("referencePoint1")
            << theParameters.value("referencePoint2")
            << ");";

        out << "                        alpha             ("
            << theParameters.value("alpha0")
            << theParameters.value("alpha1")
            << theParameters.value("alpha2")
            << ");";

        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for w component";
        out << "        Lw_Profile";
        out << "        {";

        out << "                        profile           " << profile << ";";

        out << "                        referenceValue    ("
            << theParameters.value("Lw0")
            << theParameters.value("Lw1")
            << theParameters.value("Lw2")
            << ");";

        out << "                        direction         ("
            << theParameters.value("dir0")
            << theParameters.value("dir1")
            << theParameters.value("dir2")
            << ");";

        out << "                        origin            ("
            << theParameters.value("x0")
            << theParameters.value("x1")
            << theParameters.value("x2")
            << ");";

        out << "                        referencePoint    ("
            << theParameters.value("referencePoint0")
            << theParameters.value("referencePoint1")
            << theParameters.value("referencePoint2")
            << ");";

        out << "                        alpha             ("
            << theParameters.value("alpha0")
            << theParameters.value("alpha1")
            << theParameters.value("alpha2")
            << ");";

        out << "        }";
        out << "";
        out << "";
        out << "// ************************************************************************* //";

        out.flush();
    }

    theFile.close();
}

void ExportWidget::on_btn_export_clicked()
{
    // time to export :)

    /*
     * the duplicate tree is not yet available
     * -- the code is made safe by diasabling the checkbox in this widget
     */

    // check backup style
    if (ui->duplicateTreeCheck->checkState()) {

        //
        // we need to duplicate the entire tree, then overwrite the inflow definition file
        //

        // duplicate tree:
        QString dirName = oldLocation.dirName() + ".orig";

        if (QDir(dirName).exists())
        {
            // that folder already exists! delete it!
            oldLocation.rename(dirName, dirName + ".bak");
        }
        newLocation = oldLocation;
        oldLocation.setPath(dirName);


        // write the new file
        //QDir theFile = ...;
        //this->exportInflowParameterFile(theFile);
    }
    else {

        //
        // we place new file into the existing file structure
        // but we do save one version of the existing file as
        // filename.orig before writing the new one
        //

        // save any existing file to .orig
        newLocation = oldLocation;
        newLocation.cd("constant");

        QString newFile = newLocation.absoluteFilePath("inflowProperties");
        QString origFile = newFile + ".orig";

        if (QFile(origFile).exists()) {
            qWarning() << "overwriting " << origFile;
            QFile::remove(origFile);
        }
        QFile::rename(newFile, origFile);

        qDebug() << "move" << newFile << origFile;

        // write the new file
        this->exportInflowParameterFile(newFile);
    }


}
