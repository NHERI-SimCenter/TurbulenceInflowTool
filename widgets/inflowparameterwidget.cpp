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

    /*
    ui->selectUniform->setChecked(true);

    if ( ui->selectUniform->isChecked() ) {

    }
    else if ( ui->selectExponential->isChecked() ) {

    }
    else {
        // we should not get here unless there is an issue with the UI
        qWarning() << "no model type selected.  Please report issue to SimCenter developers";
    }
    */

    ui->PHI11->setValue(0.1);
    ui->PHI21->setValue(0.0);
    ui->PHI31->setValue(0.0);
    ui->PHI22->setValue(0.1);
    ui->PHI32->setValue(0.0);
    ui->PHI33->setValue(0.1);

    ui->alpha1->setValue(0.1);
    ui->alpha2->setValue(0.1);
    ui->alpha3->setValue(0.1);

    // NEED MORE ...
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
    // collect data
    QMap<QString, double> data;
    data.clear();

    // populate data map
    double val= double(ui->modelSelectionCBX->currentIndex());
    data.insert("profile",val);

    data.insert("vel0",ui->vel->value());
    data.insert("refAngleU",ui->refAngleU->value());
    data.insert("refDistU",ui->refDistU->value());
    data.insert("alphaU",ui->alphaU->value());

    data.insert("alpha0",ui->alpha1->value());
    data.insert("alpha1",ui->alpha2->value());
    data.insert("alpha2",ui->alpha3->value());

    data.insert("phi00",ui->PHI11->value());
    data.insert("phi10",ui->PHI21->value());
    data.insert("phi20",ui->PHI31->value());
    data.insert("phi11",ui->PHI22->value());
    data.insert("phi21",ui->PHI23->value());
    data.insert("phi22",ui->PHI33->value());

    data.insert("Lu0",ui->Lux->value());
    data.insert("Lu10",ui->LuyLux->value());
    data.insert("Lu20",ui->LuzLux->value());

    data.insert("Lv0",ui->Lvx->value());
    data.insert("Lv10",ui->LvyLvx->value());
    data.insert("Lv20",ui->LvzLvx->value());

    data.insert("Lw0",ui->Lwx->value());
    data.insert("Lw10",ui->LwyLwx->value());
    data.insert("Lw20",ui->LwzLwx->value());

    data.insert("LuAlpha",ui->LuAlpha->value());
    data.insert("LvAlpha",ui->LvAlpha->value());
    data.insert("LwAlpha",ui->LwAlpha->value());

    data.insert("LuRefAngle",ui->refAngleLu->value());
    data.insert("LvRefAngle",ui->refAngleLv->value());
    data.insert("LwRefAngle",ui->refAngleLw->value());

    data.insert("LuRefDist",ui->refDistLu->value());
    data.insert("LvRefDist",ui->refDistLv->value());
    data.insert("LwRefDist",ui->refDistLw->value());

    // send the parameter map
    emit parametersReady(data);
}

