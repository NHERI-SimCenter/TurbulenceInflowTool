#include "exportwidget.h"
#include "ui_exportwidget.h"

#include <QFile>
#include <QFileInfo>
#include <QTextStream>
#include <QMap>
#include <QDir>
#include <QStandardItem>
#include <QStandardItemModel>

#include <QDebug>

ExportWidget::ExportWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::ExportWidget)
{
    ui->setupUi(this);
    theParameters.clear();
    hasParameters = false;

    UFileHead = "";
    UFileTail = "";
    clearBoundaryMap();

    ui->duplicateTreeCheck->hide();
}

ExportWidget::~ExportWidget()
{
    delete ui;
}

void ExportWidget::clearBoundaryMap(void)
{
    foreach (QString s, boundaries.keys())
    {
        if (boundaries.value(s) != nullptr) {
            delete boundaries.value(s);
        }
        boundaries.remove(s);
    }
    //qDebug() << boundaries;
}

void ExportWidget::setLocationAvailable(bool status, QDir &loc)
{
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
    case 0: { profile="uniform"; break; }
    case 1: { profile="exponential"; break; }
    default: { profile="uniform"; break; }
    }

    QFile theFile(fileName);
    if (theFile.open(QFile::WriteOnly | QFile::Truncate)) {
        QTextStream out(&theFile);

        out << "/*--------------------------------*- C++ -*----------------------------------*\\" << endl;
        out << "  =========                 |" << endl;
        out << "  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox" << endl;
        out << "   \\    /   O peration     | Website:  https://openfoam.org" << endl;
        out << "    \\  /    A nd           | Version:  6" << endl;
        out << "     \\/     M anipulation  |" << endl;
        out << "\\*---------------------------------------------------------------------------*/" << endl;
        out << "FoamFile" << endl;
        out << "{" << endl;
        out << "    version     2.0;" << endl;
        out << "    format      ascii;" << endl;
        out << "    class       dictionary;" << endl;
        out << "    location    \"constant\";" << endl;
        out << "    object      inflowProperties;" << endl;
        out << "}" << endl;
        out << "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //" << endl;
        out << "" << endl;


        out << "// mean velocity" << endl;
        out << "UDict" << endl;
        out << "{" << endl;
        out << "    referenceValue          " << theParameters.value("vel0") << ";" << endl;

        out << "    profile                 " << profile << ";" << endl;

        if ( int(theParameters.value("profile")) > 0 ) {
            out << "    referenceAngl           " << theParameters.value("refAngleU") << ";" << endl;
            out << "    referenceDist           " << theParameters.value("refDistU") << ";" << endl;
            out << "    alpha                   " << theParameters.value("alphaU") << ";" << endl;
        }
        out << "}" << endl;

        out << endl;


        out << "// turbulent intensity (symmTensorField)" << endl;
        out << "IDict" << endl;
        out << "{" << endl;
        out << "    referenceValue         ("
            << theParameters.value("phi00") << "  "
            << theParameters.value("phi10") << "  "
            << theParameters.value("phi20") << "  "
            << theParameters.value("phi11") << "  "
            << theParameters.value("phi21") << "  "
            << theParameters.value("phi22")
            << ");" << endl;

        out << "    profile                 " << profile << ";" << endl;

        if ( int(theParameters.value("profile")) > 0 ) {
            out << "    referenceAngl           " << theParameters.value("refAnglePHI") << ";" << endl;
            out << "    referenceDist           " << theParameters.value("refDistPHI") << ";" << endl;
            out << "    alpha                     ("
                << theParameters.value("alpha0") << "  "
                << theParameters.value("alpha1") << "  "
                << theParameters.value("alpha2")
                << ");" << endl;
        }

        out << "}" << endl;

        out << endl;


        out << "// turbulence length scale profile for u component" << endl;
        out << "LuxDict" << endl;
        out << "{" << endl;
        out << "    referenceValue          " << theParameters.value("Lux") << ";" << endl;

        out << "    profile                 " << profile << ";" << endl;

        if ( int(theParameters.value("profile")) > 0 ) {
            out << "    referenceAngl           " << theParameters.value("LuRefAngle") << ";" << endl;
            out << "    referenceDist           " << theParameters.value("LuRefDist") << ";" << endl;
            out << "    alpha                   " << theParameters.value("LuAlpha") << ";" << endl;
        }
        out << "}" << endl;

        out << endl;

        out << "// turbulence length scale profile for v component" << endl;
        out << "LvxDict" << endl;
        out << "{" << endl;
        out << "    referenceValue          " << theParameters.value("Luv") << ";" << endl;

        out << "    profile                 " << profile << ";" << endl;

        if ( int(theParameters.value("profile")) > 0 ) {
            out << "    referenceAngl           " << theParameters.value("LvRefAngle") << ";" << endl;
            out << "    referenceDist           " << theParameters.value("LvRefDist") << ";" << endl;
            out << "    alpha                   " << theParameters.value("LvAlpha") << ";" << endl;
        }
        out << "}" << endl;

        out << endl;


        out << "// turbulence length scale profile for w component" << endl;
        out << "LwxDict" << endl;
        out << "{" << endl;
        out << "    referenceValue          " << theParameters.value("Luw") << ";" << endl;

        out << "    profile                 " << profile << ";" << endl;

        if ( int(theParameters.value("profile")) > 0 ) {
            out << "    referenceAngl           " << theParameters.value("LwRefAngle") << ";" << endl;
            out << "    referenceDist           " << theParameters.value("LwRefDist") << ";" << endl;
            out << "    alpha                   " << theParameters.value("LwAlpha") << ";" << endl;
        }
        out << "}" << endl;

        out << endl;

        out << "LuyToLuxRatio              " << theParameters.value("Lu10") << ";" << endl;
        out << "LuzToLuxRatio              " << theParameters.value("Lu20") << ";" << endl;
        out << "LvyToLvxRatio              " << theParameters.value("Lv10") << ";" << endl;
        out << "LvzToLvxRatio              " << theParameters.value("Lv20") << ";" << endl;
        out << "LwyToLwxRatio              " << theParameters.value("Lw10") << ";" << endl;
        out << "LwzToLwxRatio              " << theParameters.value("Lw20") << ";" << endl;

        out << endl;
        out << endl;
        out << "// ************************************************************************* //" << endl;

        out.flush();
    }

    theFile.close();
}

void ExportWidget::exportUFile(QString fileName)
{
    // get the boundary condition to generate
    QString BCselected = ui->boundarySelection->currentText();

    // file handle for the U file
    QFile UFile(fileName);
    UFile.open(QFile::WriteOnly);
    QTextStream out(&UFile);

    out << UFileHead;

    foreach (QString key, boundaries.keys())
    {
        out << "    " << key << endl;
        out << "    {" << endl;

        if (key == BCselected)
        {
            QMap<QString, QString> theMap = *boundaries.value(key);

            switch (int(theParameters.value("FilterMethod"))) {
            case 0:
                out << "        type               digitalFilter;" << endl;
                switch (int(theParameters.value("shapeFunction"))) {
                case 0:
                    out << "        filterShape        gaussian;" << endl;
                    break;
                case 1:
                    out << "        filterShape        exponential;" << endl;
                    break;
                default:
                    out << "        filterShape        exponential;" << endl;
                }
                out << "        filterFactor       " << theParameters.value("filterFactor") << ";" << endl;
                out << "        gridFactor         " << theParameters.value("gridFactor") << ";" << endl;

                break;
            case 1:
                out << "        type        syntheticEddie;" << endl;
                switch (int(theParameters.value("shapeFunction"))) {
                case 0:
                    out << "        filterShape        gaussian;" << endl;
                    break;
                case 1:
                    out << "        filterShape        tent;" << endl;
                    break;
                case 2:
                    out << "        filterShape        step;" << endl;
                    break;
                default:
                    out << "        filterShape        gaussian;" << endl;
                }
                out << "        eddieDensity       " << theParameters.value("eddieDensity") << ";" << endl;

                break;
            default:
                qWarning() << "unknown turbulence model";
            }

            out << "        intersection       ( "
                << theParameters.value("intersection0") << " "
                << theParameters.value("intersection1") << " "
                << theParameters.value("intersection2") << " );" << endl;
            out << "        yOffset            " << theParameters.value("yOffset") << ";" << endl;
            out << "        zOffset            " << theParameters.value("zOffset") << ";" << endl;

            if (theMap.contains("type"))         theMap.remove("type");
            if (theMap.contains("filterShape"))  theMap.remove("filterShape");
            if (theMap.contains("filterFactor")) theMap.remove("filterFactor");
            if (theMap.contains("gridFactor"))   theMap.remove("gridFactor");
            if (theMap.contains("eddieDensity")) theMap.remove("eddieDensity");

            if (theMap.contains("intersection"))    theMap.remove("intersection");
            if (theMap.contains("yOffset"))         theMap.remove("yOffset");
            if (theMap.contains("zOffset"))         theMap.remove("zOffset");

            foreach (QString s, theMap.keys() )
            {
                out << "        " << s << "    " << theMap.value(s) << ";" << endl;
            }
        }
        else {
            foreach (QString s, (boundaries.value(key))->keys() )
            {
                out << "        " << s << "    " << (boundaries.value(key))->value(s) << ";" << endl;
            }
        }
        out << "    }" << endl;
        out << endl;
    }

    out << UFileTail;

    UFile.close();
}

void ExportWidget::exportControlDictFile(QString fileName)
{
    // file handle for the controlDict file
    QFile CDict(fileName);
    CDict.open(QFile::ReadOnly);
    CDictContents = CDict.readAll();
    CDict.close();

    CDict.open(QFile::WriteOnly);
    QTextStream out(&CDict);

    QList<QByteArray> CDictList = CDictContents.split('\n');
    foreach (QByteArray line, CDictList)
    {
        if (line.contains("")) {
            out << "libs" << endl;
            out << "(" << endl;
            out << "    \"libturbulentInflow.so\"" << endl;
            out << ");" << endl;
            out << endl;
        }

        out << line << endl;
    }

    CDict.close();
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

        // we place new file into the existing file structure
        // but we do save one version of the existing file as
        // filename.orig before writing the new one

        //
        // ... inflowProperties file
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

        //
        // ... U file
        //

        newLocation = oldLocation;
        newLocation.cd("0");

        newFile  = newLocation.absoluteFilePath("U");
        origFile = newFile + ".orig";

        if (QFile(origFile).exists()) {
            qWarning() << "overwriting " << origFile;
            QFile::remove(origFile);
        }
        QFile::rename(newFile, origFile);

        qDebug() << "move" << newFile << origFile;

        // update U file
        this->exportUFile(newFile);

        //
        // ... controlDict file
        //

        newLocation = oldLocation;
        newLocation.cd("systen");

        newFile  = newLocation.absoluteFilePath("controlDict");
        origFile = newFile + ".orig";

        if (QFile(origFile).exists()) {
            qWarning() << "overwriting " << origFile;
            QFile::remove(origFile);
        }
        QFile::rename(newFile, origFile);

        qDebug() << "move" << newFile << origFile;

        // update controlDict file
        this->exportControlDictFile(newFile);
    }
}

void ExportWidget::setUFileData(QByteArray &head, QByteArray &tail, QMap<QString, QMap<QString, QString> * > &data)
{
    UFileHead = head;
    UFileTail = tail;
    boundaries = data;
}

void ExportWidget::setBoundarySelection(int index)
{
    ui->boundarySelection->setCurrentIndex(index);
}

void ExportWidget::on_boundarySelection_currentIndexChanged(int index)
{
    emit boundarySelection(index);
}

void ExportWidget::setModel(QStandardItemModel *theModel)
{
    ui->boundarySelection->setModel(theModel);
    validSourcePresent = true;
}
