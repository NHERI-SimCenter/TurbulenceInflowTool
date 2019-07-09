#include "inflowparameterwidget.h"
#include "ui_inflowparameterwidget.h"

#include "math.h"

#include <QFileDialog>
#include <QDebug>

InflowParameterWidget::InflowParameterWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::InflowParameterWidget)
{
    ui->setupUi(this);
    ui->sourceGroup->hide();
    setDefaultParameters();
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
    theParameters["refDistU"] = 0.0;
    theParameters["alphaU"] = 0.0;

    theParameters["alpha0"] = 0.0;
    theParameters["alpha1"] = 0.0;
    theParameters["alpha2"] = 0.0;

    theParameters["phi00"] = 0.1;
    theParameters["phi10"] = 0.0;
    theParameters["phi20"] = 0.0;
    theParameters["phi11"] = 0.1;
    theParameters["phi21"] = 0.0;
    theParameters["phi22"] = 0.1;

    theParameters["Lu0"] = 1.0;
    theParameters["Lu10"] = 0.0;
    theParameters["Lu20"] = 0.0;

    theParameters["Lv0"] = 1.0;
    theParameters["Lv10"] = 0.0;
    theParameters["Lv20"] = 0.0;

    theParameters["Lw0"] = 1.0;
    theParameters["Lw10"] = 0.0;
    theParameters["Lw20"] = 0.0;

    theParameters["LuAlpha"] = 0.0;
    theParameters["LvAlpha"] = 0.0;
    theParameters["LwAlpha"] = 0.0;

    theParameters["LuRefAngle"] = 0.0;
    theParameters["LvRefAngle"] = 0.0;
    theParameters["LwRefAngle"] = 0.0;

    theParameters["LuRefDist"] = 0.0;
    theParameters["LvRefDist"] = 0.0;
    theParameters["LwRefDist"] = 0.0;

    /* for use in U file */

    theParameters["FilterMethod"] = 0;

    theParameters["shapeFunction"] = 0;
    theParameters["gridFactor"] = 1.0;
    theParameters["filterFactor"] = 4;

    theParameters["velocityShape"] = 0;
    theParameters["eddyDensity"] = 0.0;

    theParameters["intersection0"] = 0.0;
    theParameters["intersection1"] = 1.0;
    theParameters["intersection2"] = 0.0;
    theParameters["yOffset"] = 0.0;
    theParameters["zOffset"] = 0.0;

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

    theParameters.insert("alpha0",ui->alpha1->value());
    theParameters.insert("alpha1",ui->alpha2->value());
    theParameters.insert("alpha2",ui->alpha3->value());

    theParameters.insert("phi00",ui->PHI11->value());
    theParameters.insert("phi10",ui->PHI21->value());
    theParameters.insert("phi20",ui->PHI31->value());
    theParameters.insert("phi11",ui->PHI22->value());
    theParameters.insert("phi21",ui->PHI23->value());
    theParameters.insert("phi22",ui->PHI33->value());

    theParameters.insert("Lu0",ui->Lux->value());
    theParameters.insert("Lu10",ui->LuyLux->value());
    theParameters.insert("Lu20",ui->LuzLux->value());

    theParameters.insert("Lv0",ui->Lvx->value());
    theParameters.insert("Lv10",ui->LvyLvx->value());
    theParameters.insert("Lv20",ui->LvzLvx->value());

    theParameters.insert("Lw0",ui->Lwx->value());
    theParameters.insert("Lw10",ui->LwyLwx->value());
    theParameters.insert("Lw20",ui->LwzLwx->value());

    theParameters.insert("LuAlpha",ui->LuAlpha->value());
    theParameters.insert("LvAlpha",ui->LvAlpha->value());
    theParameters.insert("LwAlpha",ui->LwAlpha->value());

    theParameters.insert("LuRefAngle",ui->refAngleLu->value());
    theParameters.insert("LvRefAngle",ui->refAngleLv->value());
    theParameters.insert("LwRefAngle",ui->refAngleLw->value());

    theParameters.insert("LuRefDist",ui->refDistLu->value());
    theParameters.insert("LvRefDist",ui->refDistLv->value());
    theParameters.insert("LwRefDist",ui->refDistLw->value());

    /* for use in U file */

    if (ui->RB_digitalFilter->isChecked())
        { theParameters.insert("FilterMethod",0); }
    else
        { theParameters.insert("FilterMethod",1); }

    theParameters.insert("shapeFunction",ui->shapeFunction->currentIndex());
    theParameters.insert("gridFactor",ui->gridFactor->value());
    theParameters.insert("filterFactor",ui->filterFactor->value());

    theParameters.insert("velocityShape",ui->velocityShape->currentIndex());
    theParameters.insert("eddyDensity",ui->eddyDensity->value());

    theParameters.insert("intersection0",ui->dir1->value());
    theParameters.insert("intersection1",ui->dir2->value());
    theParameters.insert("intersection2",ui->dir3->value());
    theParameters.insert("yOffset",ui->yOffset->value());
    theParameters.insert("zOffset",ui->zOffset->value());

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

    ui->alpha1->setValue(theParameters.value("alpha0"));
    ui->alpha2->setValue(theParameters.value("alpha1"));
    ui->alpha3->setValue(theParameters.value("alpha2"));

    ui->PHI11->setValue(theParameters.value("phi00"));
    ui->PHI21->setValue(theParameters.value("phi10"));
    ui->PHI31->setValue(theParameters.value("phi20"));
    ui->PHI22->setValue(theParameters.value("phi11"));
    ui->PHI23->setValue(theParameters.value("phi21"));
    ui->PHI33->setValue(theParameters.value("phi22"));

    ui->Lux->setValue(theParameters.value("Lu0"));
    ui->LuyLux->setValue(theParameters.value("Lu10"));
    ui->LuzLux->setValue(theParameters.value("Lu20"));

    ui->Lvx->setValue(theParameters.value("Lv0"));
    ui->LvyLvx->setValue(theParameters.value("Lv10"));
    ui->LvzLvx->setValue(theParameters.value("Lv20"));

    ui->Lwx->setValue(theParameters.value("Lw0"));
    ui->LwyLwx->setValue(theParameters.value("Lw10"));
    ui->LwzLwx->setValue(theParameters.value("Lw20"));

    ui->LuAlpha->setValue(theParameters.value("LuAlpha"));
    ui->LvAlpha->setValue(theParameters.value("LvAlpha"));
    ui->LwAlpha->setValue(theParameters.value("LwAlpha"));

    ui->refAngleLu->setValue(theParameters.value("LuRefAngle"));
    ui->refAngleLv->setValue(theParameters.value("LvRefAngle"));
    ui->refAngleLw->setValue(theParameters.value("LwRefAngle"));

    ui->refDistLu->setValue(theParameters.value("LuRefDist"));
    ui->refDistLv->setValue(theParameters.value("LvRefDist"));
    ui->refDistLw->setValue(theParameters.value("LwRefDist"));

    /* for use in U file */

    ui->RB_digitalFilter->setChecked(int(theParameters.value("FilterMethod"))==0?true:false);

    ui->shapeFunction->setCurrentIndex(int(theParameters.value("shapeFunction")));
    ui->gridFactor->setValue(theParameters.value("gridFactor"));
    ui->filterFactor->setValue(int(theParameters.value("filterFactor")));

    ui->velocityShape->setCurrentIndex(int(theParameters.value("velocityShape")));
    ui->eddyDensity->setValue(theParameters.value("eddyDensity"));

    ui->dir1->setValue(theParameters.value("intersection0"));
    ui->dir2->setValue(theParameters.value("intersection1"));
    ui->dir3->setValue(theParameters.value("intersection2"));
    ui->yOffset->setValue(theParameters.value("yOffset"));
    ui->zOffset->setValue(theParameters.value("zOffset"));
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
    ui->lengthScaleGroup->show();

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

    ui->label_referenceDist->hide();
    ui->label_referenceAngle->hide();
    ui->label_alphas->hide();

    ui->refDistLu->hide();
    ui->refAngleLu->hide();
    ui->LuAlpha->hide();
    ui->refDistLv->hide();
    ui->refAngleLv->hide();
    ui->LvAlpha->hide();
    ui->refDistLw->hide();
    ui->refAngleLw->hide();
    ui->LwAlpha->hide();

    ui->line3->hide();
    ui->line4->hide();
}

void InflowParameterWidget::setExponentialTurbulent(void)
{
    ui->velocityGroup->show();
    ui->phiTensorGroup->show();
    ui->lengthScaleGroup->show();

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

    ui->label_referenceDist->show();
    ui->label_referenceAngle->show();
    ui->label_alphas->show();

    ui->refDistLu->show();
    ui->refAngleLu->show();
    ui->LuAlpha->show();
    ui->refDistLv->show();
    ui->refAngleLv->show();
    ui->LvAlpha->show();
    ui->refDistLw->show();
    ui->refAngleLw->show();
    ui->LwAlpha->show();

    ui->line3->show();
    ui->line4->show();
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
    ui->stackedMethods->setCurrentIndex((ui->RB_digitalFilter->isChecked())?0:1);
}

void InflowParameterWidget::on_RB_syntheticEddie_clicked()
{
    ui->stackedMethods->setCurrentIndex((ui->RB_digitalFilter->isChecked())?0:1);
}
