#include "exportwidget.h"
#include "ui_exportwidget.h"

#include <QFile>
#include <QTextStream>
#include <QMap>

#include <QDebug>

ExportWidget::ExportWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::ExportWidget)
{
    ui->setupUi(this);
    theParameters = nullptr;
}

ExportWidget::~ExportWidget()
{
    delete ui;
}

void ExportWidget::setLocationAvailable(bool status, QDir &loc)
{
    qDebug() << "exportWidget: signal received";

    if (status) {
        hasLocation = false;
        oldLocation = QDir(".");
        newLocation = QDir(".");
    }
    else {
        hasLocation = true;
        oldLocation = loc;
        newLocation = QDir(".");
    }
}

void ExportWidget::exportInflowParameterFile(QString fileName)
{
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
        out << "                        referenceValue    (10 0 0);";
        out << "                        profile           exponential;";
        out << "                        direction         (0 0 1);";
        out << "                        origin            (0 0 0);";

        out << "                        referencePoint    (";
        out << 0 << 0 << 0.364 ;
        out << ");";

        out << "                        alpha             (0.326 0 0);";
        out << "        }";
        out << "";
        out << "";
        out << "        // turbulent intensity (symmTensorField)";
        out << "        I_Profile";
        out << "        {";
        out << "                        profile           exponential;";
        out << "                        referenceValue    (0.208 0 0 0.182 0 0.152);";
        out << "                        direction         (0 0 1);";
        out << "                        origin            (0 0 0);                 ";
        out << "                        referencePoint    (0 0 0.364);";
        out << "                        alpha             (-0.191 0 0 -0.123 0 -0.005);";
        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for u component";
        out << "        Lu_Profile";
        out << "        {";
        out << "                        profile           exponential;";
        out << "                        referenceValue    (0.302 0.302 0.302);";
        out << "                        direction         (0 0 1);";
        out << "                        origin            (0 0 0);                 ";
        out << "                        referencePoint    (0 0 0.254);      ";
        out << "                        alpha             (0.473 0.473 0.473);         ";
        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for v component";
        out << "        Lv_Profile";
        out << "        {";
        out << "                        profile           exponential;";
        out << "                        referenceValue    (0.0815 0.0815 0.0815);";
        out << "                        direction         (0 0 1);";
        out << "                        origin            (0 0 0);                 ";
        out << "                        referencePoint    (0 0 0.254);";
        out << "                        alpha             (0.881 0.881 0.881);";
        out << "        }";
        out << "";
        out << "        // turbulence length scale profile for w component";
        out << "        Lw_Profile";
        out << "        {";
        out << "                        profile           exponential;";
        out << "                        referenceValue    (0.0326 0.0326 0.0326);";
        out << "                        direction         (0 0 1);";
        out << "                        origin            (0 0 0);                 ";
        out << "                        referencePoint    (0 0 0.254);";
        out << "                        alpha             (1.539 1.539 1.539);";
        out << "        }";
        out << "";
        out << "";
        out << "// ************************************************************************* //";

        out.flush();
    }

    theFile.close();
}
