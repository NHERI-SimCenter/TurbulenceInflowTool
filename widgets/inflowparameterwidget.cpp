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

    ui->z01->setValue(0.0);
    ui->z02->setValue(0.0);
    ui->z03->setValue(0.0);

    ui->zVector1->setValue(0.0);
    ui->zVector1->setValue(0.0);
    ui->zVector1->setValue(1.0);

    ui->alpha1->setValue(0.1);
    ui->alpha2->setValue(0.1);
    ui->alpha3->setValue(0.1);
}

void InflowParameterWidget::on_btnNormalize_clicked()
{
    double z1 = ui->zVector1->value();
    double z2 = ui->zVector2->value();
    double z3 = ui->zVector3->value();

    double norm = sqrt(z1*z1 + z2*z2 + z3+z3);

    if (norm > 1.0e-3)
    {
        ui->zVector1->setValue(z1/norm);
        ui->zVector2->setValue(z2/norm);
        ui->zVector3->setValue(z3/norm);
    }
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
        this->setLinearLaminar();
        break;
    case 1:
        this->setExponentialLaminar();
        break;
    case 2:
        this->setLinearTurbulent();
        break;
    case 3:
        this->setExponentialTurbulent();
        break;
    default:
        qWarning() << "Unknown boundary condition type selected" ;
    }
}

void InflowParameterWidget::setLinearLaminar(void)
{
    ui->phiTensorGroup->hide();
    ui->alphaParameterGroup->hide();
    ui->lengthScaleGroup->hide();
    ui->referencePointGroup->hide();
    ui->velocityGroup->show();
}

void InflowParameterWidget::setExponentialLaminar(void)
{
    ui->phiTensorGroup->show();
    ui->alphaParameterGroup->show();
    ui->lengthScaleGroup->hide();
    ui->referencePointGroup->show();
    ui->velocityGroup->hide();
}

void InflowParameterWidget::setLinearTurbulent(void)
{
    ui->phiTensorGroup->hide();
    ui->alphaParameterGroup->hide();
    ui->lengthScaleGroup->show();
    ui->referencePointGroup->show();
    ui->velocityGroup->show();
}

void InflowParameterWidget::setExponentialTurbulent(void)
{
    ui->phiTensorGroup->show();
    ui->alphaParameterGroup->show();
    ui->lengthScaleGroup->show();
    ui->referencePointGroup->show();
    ui->velocityGroup->hide();
}

void InflowParameterWidget::sendParameterMap(void)
{
    // collect data
    QMap<QString, double> data;
    data.clear();

    // populate data map
    double val= double(ui->modelSelectionCBX->currentIndex());
    data.insert("profile",val);

    data.insert("x0",ui->z01->value());
    data.insert("x1",ui->z02->value());
    data.insert("x2",ui->z03->value());

    data.insert("dir0",ui->zVector1->value());
    data.insert("dir1",ui->zVector2->value());
    data.insert("dir2",ui->zVector3->value());

    data.insert("vel0",ui->vel1->value());
    data.insert("vel1",ui->vel2->value());
    data.insert("vel2",ui->vel3->value());

    data.insert("alpha0",ui->alpha1->value());
    data.insert("alpha1",ui->alpha2->value());
    data.insert("alpha2",ui->alpha3->value());

    data.insert("phi00",ui->PHI11->value());
    data.insert("phi10",ui->PHI21->value());
    data.insert("phi20",ui->PHI31->value());
    data.insert("phi11",ui->PHI22->value());
    data.insert("phi21",ui->PHI23->value());
    data.insert("phi22",ui->PHI33->value());

    data.insert("referencePoint0",ui->tVector1->value());
    data.insert("referencePoint1",ui->tVector2->value());
    data.insert("referencePoint2",ui->tVector3->value());

    data.insert("Lu0",ui->Lux->value());
    data.insert("Lu1",ui->Luy->value());
    data.insert("Lu2",ui->Luz->value());

    data.insert("Lv0",ui->Lvx->value());
    data.insert("Lv1",ui->Lvy->value());
    data.insert("Lv2",ui->Lvz->value());

    data.insert("Lw0",ui->Lwx->value());
    data.insert("Lw1",ui->Lwy->value());
    data.insert("Lw2",ui->Lwz->value());

    // send the parameter map
    emit parametersReady(data);
}

