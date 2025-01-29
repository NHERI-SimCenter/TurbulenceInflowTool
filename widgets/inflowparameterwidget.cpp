#include "inflowparameterwidget.h"
#include "ui_inflowparameterwidget.h"

#include "math.h"

#include <QFileDialog>
#include <QJsonObject>
#include <QDebug>

InflowParameterWidget::InflowParameterWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::InflowParameterWidget)
{
    ui->setupUi(this);
    ui->sourceGroup->hide();
    setDefaultParameters();

    theParameters.clear();
    hasParameters = false;
}

InflowParameterWidget::~InflowParameterWidget()
{
    delete ui;
}

void InflowParameterWidget::selectSourceLocation(void)
{
    QFileDialog *dlg = new QFileDialog();
    dlg->setFileMode(QFileDialog::Directory);
    if (dlg->exec())
    {
        QDir sourceFolder = dlg->directory();
        ui->sourceLocationDisplay->setText(sourceFolder.canonicalPath());

        qDebug() << sourceFolder;
    }
    delete dlg;

}

void InflowParameterWidget::setDefaultParameters()
{
    this->on_modelSelectionCBX_currentIndexChanged(0);

    theParameters.clear();

    /* for use in inflowProperties file */

    theParameters["profile"] = 0;

    theParameters["vel0"] = 1.0;
    theParameters["refAngleU"] = 0.0;
    theParameters["refDistU"] = 1.0;
    theParameters["alphaU"] = 0.0;

    theParameters["refAnglePHI"] = 0.0;
    theParameters["refDistPHI"] = 1.0;

    theParameters["alpha0"] = 0.0;
    theParameters["alpha1"] = 0.0;
    theParameters["alpha2"] = 0.0;

    theParameters["phi00"] = 1.0;
    theParameters["phi10"] = 0.0;
    theParameters["phi20"] = 0.0;
    theParameters["phi11"] = 1.0;
    theParameters["phi21"] = 0.0;
    theParameters["phi22"] = 1.0;
    
    theParameters["refAngleL"] = 0.0;
    theParameters["refDistL"] = 1.0;
    
    theParameters["L11"] = 1.0;
    theParameters["L12"] = 1.0;
    theParameters["L13"] = 1.0;
    theParameters["L21"] = 1.0;
    theParameters["L22"] = 1.0;
    theParameters["L23"] = 1.0;
    theParameters["L31"] = 1.0;
    theParameters["L32"] = 1.0;
    theParameters["L33"] = 1.0;
    
    theParameters["alpha11"] = 0.0;
    theParameters["alpha12"] = 0.0;
    theParameters["alpha13"] = 0.0;
    theParameters["alpha21"] = 0.0;
    theParameters["alpha22"] = 0.0;
    theParameters["alpha23"] = 0.0;
    theParameters["alpha31"] = 0.0;
    theParameters["alpha32"] = 0.0;
    theParameters["alpha33"] = 0.0;

    /* for use in U file */

    theParameters["FilterMethod"] = 0;

    theParameters["filterType"] = 0;
    theParameters["gridFactor"] = 1.0;
    theParameters["filterFactor"] = 2;

    theParameters["eddyType"] = 0;
    theParameters["eddyDensity"] = 0.0;

    theParameters["intersection0"] = 0.0;
    theParameters["intersection1"] = 0.0;
    theParameters["intersection2"] = 0.0;
    
    theParameters["offset0"] = 0.0;
    theParameters["offset1"] = 0.0;
    theParameters["offset2"] = 0.0;

    hasParameters = true;

    refreshDisplay();
}

void InflowParameterWidget::refreshParameterMap(void)
{
    // collect data
    theParameters.clear();

    //
    // populate theParameters map
    //

    /* for use in inflowProperties file */

    theParameters.insert("profile",double(ui->modelSelectionCBX->currentIndex()));

    theParameters.insert("vel0",ui->vel->value());
    theParameters.insert("refAngleU",ui->refAngleU->value());
    theParameters.insert("refDistU",ui->refDistU->value());
    theParameters.insert("alphaU",ui->alphaU->value());

    theParameters.insert("refAnglePHI",ui->refAnglePHI->value());
    theParameters.insert("refDistPHI",ui->refDistPHI->value());

    theParameters.insert("alpha0",ui->alpha1->value());
    theParameters.insert("alpha1",ui->alpha2->value());
    theParameters.insert("alpha2",ui->alpha3->value());

    theParameters.insert("phi00",ui->PHI11->value());
    theParameters.insert("phi10",ui->PHI21->value());
    theParameters.insert("phi20",ui->PHI31->value());
    theParameters.insert("phi11",ui->PHI22->value());
    theParameters.insert("phi21",ui->PHI23->value());
    theParameters.insert("phi22",ui->PHI33->value());
    
    theParameters.insert("refAngleL",ui->refAngleL->value());
    theParameters.insert("refDistL",ui->refDistL->value());

    theParameters.insert("alpha11",ui->alpha11->value());
    theParameters.insert("alpha12",ui->alpha12->value());
    theParameters.insert("alpha13",ui->alpha13->value());
    theParameters.insert("alpha21",ui->alpha21->value());
    theParameters.insert("alpha22",ui->alpha22->value());
    theParameters.insert("alpha23",ui->alpha23->value());
    theParameters.insert("alpha31",ui->alpha31->value());
    theParameters.insert("alpha32",ui->alpha32->value());
    theParameters.insert("alpha33",ui->alpha33->value());

    theParameters.insert("L11",ui->L11->value());
    theParameters.insert("L12",ui->L12->value());
    theParameters.insert("L13",ui->L13->value());
    theParameters.insert("L21",ui->L21->value());
    theParameters.insert("L22",ui->L22->value());
    theParameters.insert("L23",ui->L23->value());
    theParameters.insert("L31",ui->L31->value());
    theParameters.insert("L32",ui->L32->value());
    theParameters.insert("L33",ui->L33->value());

    /* for use in U file */

    // there must be four options FIX IT!

    if (ui->RB_digitalFilter->isChecked())
        { theParameters.insert("FilterMethod",0); }
    else if (ui->RB_syntheticEddie->isChecked())
        { theParameters.insert("FilterMethod",1); }
    else if (ui->RB_divergenceFree->isChecked())
        { theParameters.insert("FilterMethod",2); }
    else if (ui->RB_turbulentSpot->isChecked())
        { theParameters.insert("FilterMethod",3); }
    else
        { theParameters.insert("FilterMethod",0); }

    theParameters.insert("filterType",ui->filterType->currentIndex());
    theParameters.insert("gridFactor",ui->gridFactor->value());
    theParameters.insert("filterFactor",ui->filterFactor->value());

    theParameters.insert("eddyType",ui->eddyType->currentIndex());
    theParameters.insert("eddyDensity",ui->eddyDensity->value());
    theParameters.insert("divergenceFreeEddyDensity",ui->divEddyDensity->value());
    theParameters.insert("turbulentSpotDensity",ui->turbulentSpotDensity->value());
    
    if (ui->RB_turbulentSpotTypeL->isChecked()) {
        theParameters.insert("turbulentSpotType", -1.0);
    }
    else {
        theParameters.insert("turbulentSpotType", 1.0);
    }

    theParameters.insert("periodicY",ui->CBx_periodicY->isChecked()?1:0);
    theParameters.insert("periodicZ",ui->CBx_periodicZ->isChecked()?1:0);
    theParameters.insert("cleanRestart",ui->CBx_cleanRestart->isChecked()?1:0);
    theParameters.insert("interpolateParameters",ui->CBx_interpolateParameters->isChecked()?1:0);

    theParameters.insert("intersection0",ui->dir1->value());
    theParameters.insert("intersection1",ui->dir2->value());
    theParameters.insert("intersection2",ui->dir3->value());
    
    theParameters.insert("offset0",ui->offset0->value());
    theParameters.insert("offset1",ui->offset1->value());
    theParameters.insert("offset2",ui->offset2->value());

    hasParameters = true;
}

void InflowParameterWidget::refreshDisplay(void)
{
    /* for use in inflowProperties file */

    ui->modelSelectionCBX->setCurrentIndex(int(theParameters.value("profile")));

    ui->vel->setValue(theParameters.value("vel0"));
    ui->refAngleU->setValue(theParameters.value("refAngleU"));
    ui->refDistU->setValue(theParameters.value("refDistU"));
    ui->alphaU->setValue(theParameters.value("alphaU"));

    ui->refAnglePHI->setValue(theParameters.value("refAnglePHI"));
    ui->refDistPHI->setValue(theParameters.value("refDistPHI"));

    ui->alpha1->setValue(theParameters.value("alpha0"));
    ui->alpha2->setValue(theParameters.value("alpha1"));
    ui->alpha3->setValue(theParameters.value("alpha2"));

    ui->PHI11->setValue(theParameters.value("phi00"));
    ui->PHI21->setValue(theParameters.value("phi10"));
    ui->PHI31->setValue(theParameters.value("phi20"));
    ui->PHI22->setValue(theParameters.value("phi11"));
    ui->PHI23->setValue(theParameters.value("phi21"));
    ui->PHI33->setValue(theParameters.value("phi22"));
    
    ui->refAngleL->setValue(theParameters.value("refAngleL"));
    ui->refDistL->setValue(theParameters.value("refDistL"));

    ui->alpha11->setValue(theParameters.value("alpha11"));
    ui->alpha12->setValue(theParameters.value("alpha12"));
    ui->alpha13->setValue(theParameters.value("alpha13"));
    ui->alpha21->setValue(theParameters.value("alpha21"));
    ui->alpha22->setValue(theParameters.value("alpha22"));
    ui->alpha23->setValue(theParameters.value("alpha23"));
    ui->alpha31->setValue(theParameters.value("alpha31"));
    ui->alpha32->setValue(theParameters.value("alpha32"));
    ui->alpha33->setValue(theParameters.value("alpha33"));
    
    ui->L11->setValue(theParameters.value("L11"));
    ui->L12->setValue(theParameters.value("L12"));
    ui->L13->setValue(theParameters.value("L13"));
    ui->L21->setValue(theParameters.value("L21"));
    ui->L22->setValue(theParameters.value("L22"));
    ui->L23->setValue(theParameters.value("L23"));
    ui->L31->setValue(theParameters.value("L31"));
    ui->L32->setValue(theParameters.value("L32"));
    ui->L33->setValue(theParameters.value("L33"));

    /* for use in U file */

    ui->RB_digitalFilter->setChecked(int(theParameters.value("FilterMethod"))==0?true:false);
    ui->RB_syntheticEddie->setChecked(int(theParameters.value("FilterMethod"))==1?true:false);
    ui->RB_divergenceFree->setChecked(int(theParameters.value("FilterMethod"))==2?true:false);
    ui->RB_turbulentSpot->setChecked(int(theParameters.value("FilterMethod"))==3?true:false);

    ui->filterType->setCurrentIndex(int(theParameters.value("filterType")));
    ui->gridFactor->setValue(theParameters.value("gridFactor"));
    ui->filterFactor->setValue(int(theParameters.value("filterFactor")));

    ui->eddyType->setCurrentIndex(int(theParameters.value("eddyType")));
    ui->eddyDensity->setValue(theParameters.value("eddyDensity"));

    ui->divEddyDensity->setValue(theParameters.value("divergenceFreeEddyDensity"));
    ui->turbulentSpotDensity->setValue(theParameters.value("turbulentSpotDensity"));
    
    if (theParameters.value("turbulentSpotType") > 0.0) {
        ui->RB_turbulentSpotTypeL->setChecked(false);
        ui->RB_turbulentSpotTypeR->setChecked(true);
    }
    else {
        ui->RB_turbulentSpotTypeL->setChecked(true);
        ui->RB_turbulentSpotTypeR->setChecked(false);
    }

    ui->CBx_periodicY->setChecked(theParameters.value("periodicY") > 0.1);
    ui->CBx_periodicZ->setChecked(theParameters.value("periodicZ") > 0.1);
    ui->CBx_cleanRestart->setChecked(theParameters.value("cleanRestart") > 0.1);
    ui->CBx_interpolateParameters->setChecked(theParameters.value("interpolateParameters") > 0.1);

    ui->dir1->setValue(theParameters.value("intersection0"));
    ui->dir2->setValue(theParameters.value("intersection1"));
    ui->dir3->setValue(theParameters.value("intersection2"));
    
    ui->offset0->setValue(theParameters.value("offset0"));
    ui->offset1->setValue(theParameters.value("offset1"));
    ui->offset2->setValue(theParameters.value("offset2"));
}

void InflowParameterWidget::on_PHI21_valueChanged(double arg1)
{
    ui->PHI12->setValue(arg1);
}

void InflowParameterWidget::on_PHI31_valueChanged(double arg1)
{
    ui->PHI13->setValue(arg1);
}

void InflowParameterWidget::on_PHI32_valueChanged(double arg1)
{
    ui->PHI23->setValue(arg1);
}

void InflowParameterWidget::on_resetButton_clicked()
{
    // set UI to default parameter values
    setDefaultParameters();
}

void InflowParameterWidget::on_modelSelectionCBX_currentIndexChanged(int index)
{
    // this is where we get a mode

    switch (index) {
    case 0:
        this->setUniformTurbulent();
        break;
    case 1:
        this->setExponentialTurbulent();
        break;
    default:
        qWarning() << "Unknown boundary condition type selected" ;
    }
}

void InflowParameterWidget::setUniformTurbulent(void)
{
    ui->velocityGroup->show();
    ui->phiTensorGroup->show();
    ui->LTensorGroup->show();

    // hide extension parameters
    ui->label_refDistU->hide();
    ui->label_refAngleU->hide();
    ui->label_alphaU->hide();

    ui->refDistU->hide();
    ui->refAngleU->hide();
    ui->alphaU->hide();

    ui->label_refDistPHI->hide();
    ui->label_refAnglePHI->hide();
    ui->label_alphaPHI->hide();

    ui->refDistPHI->hide();
    ui->refAnglePHI->hide();
    ui->alpha1->hide();
    ui->alpha2->hide();
    ui->alpha3->hide();
    
    ui->label_refDistL->hide();
    ui->label_refAngleL->hide();
    ui->label_alphaL->hide();

    ui->refDistL->hide();
    ui->refAngleL->hide();
    ui->alpha11->hide();
    ui->alpha12->hide();
    ui->alpha13->hide();
    ui->alpha21->hide();
    ui->alpha22->hide();
    ui->alpha23->hide();
    ui->alpha31->hide();
    ui->alpha32->hide();
    ui->alpha33->hide();
}

void InflowParameterWidget::setExponentialTurbulent(void)
{
    ui->velocityGroup->show();
    ui->phiTensorGroup->show();
    ui->LTensorGroup->show();

    // show extension parameters
    ui->label_refDistU->show();
    ui->label_refAngleU->show();
    ui->label_alphaU->show();

    ui->refDistU->show();
    ui->refAngleU->show();
    ui->alphaU->show();

    ui->label_refDistPHI->show();
    ui->label_refAnglePHI->show();
    ui->label_alphaPHI->show();

    ui->refDistPHI->show();
    ui->refAnglePHI->show();
    ui->alpha1->show();
    ui->alpha2->show();
    ui->alpha3->show();
    
    ui->label_refDistL->show();
    ui->label_refAngleL->show();
    ui->label_alphaL->show();

    ui->refDistL->show();
    ui->refAngleL->show();
    ui->alpha11->show();
    ui->alpha12->show();
    ui->alpha13->show();
    ui->alpha21->show();
    ui->alpha22->show();
    ui->alpha23->show();
    ui->alpha31->show();
    ui->alpha32->show();
    ui->alpha33->show();
}

void InflowParameterWidget::sendParameterMap(void)
{
    this->refreshParameterMap();

    // collect data
    QMap<QString, double> data(theParameters);

    // send the parameter map
    emit parametersReady(data);
}


void InflowParameterWidget::on_RB_digitalFilter_clicked()
{
    ui->stackedMethods->setCurrentIndex(0);
}

void InflowParameterWidget::on_RB_syntheticEddie_clicked()
{
    ui->stackedMethods->setCurrentIndex(1);
}

void InflowParameterWidget::on_RB_divergenceFree_clicked()
{
    ui->stackedMethods->setCurrentIndex(2);
}

void InflowParameterWidget::on_RB_turbulentSpot_clicked()
{
    ui->stackedMethods->setCurrentIndex(3);
}


bool InflowParameterWidget::outputToJSON(QJsonObject &rvObject)
{
    refreshParameterMap();

    // just need to send the class type here.. type needed in object in case user screws up
    rvObject["type"]="CFD-Inflow";

    rvObject["EventClassification"]="Wind";

    foreach (QString key, theParameters.keys())
    {
        rvObject[key] = theParameters.value(key);
    }

    return true;
}

bool InflowParameterWidget::inputFromJSON(QJsonObject &rvObject)
{
    // initialize theParameters to reflect all properties
    refreshParameterMap();

    // update theParameters using information from the JSON file
    foreach (QString key, theParameters.keys())
    {
        if (rvObject.contains(key)) {
          QJsonValue theValue = rvObject[key];
          theParameters[key] = theValue.toDouble();
        }
        else
          return false;
    }

    // update parameter values
    refreshDisplay();

    return true;
}

void InflowParameterWidget::reset(void)
{
    setDefaultParameters();
}

void InflowParameterWidget::on_CBx_interpolateParameters_clicked()
{
    bool interpolate = ui->CBx_interpolateParameters->isChecked() ;

    ui->localCoordinateSystemGroup->setEnabled(!interpolate);
    ui->typeSelectionGroup->setEnabled(!interpolate);
    ui->velocityGroup->setEnabled(!interpolate);
    ui->phiTensorGroup->setEnabled(!interpolate);
    ui->LTensorGroup->setEnabled(!interpolate);
}
