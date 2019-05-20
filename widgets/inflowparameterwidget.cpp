#include "inflowparameterwidget.h"
#include "ui_inflowparameterwidget.h"

#include "math.h"

#include <QDebug>

InflowParameterWidget::InflowParameterWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::InflowParameterWidget)
{
    ui->setupUi(this);
    setDefaultParameters();
}

InflowParameterWidget::~InflowParameterWidget()
{
    delete ui;
}

void InflowParameterWidget::setDefaultParameters()
{
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

void InflowParameterWidget::on_selectUniform_clicked(bool checked)
{
    // adjust UI for uniform distribution
    if ( checked ) {

        // deactivate alpha parameters

        ui->alpha1->setEnabled(false);
        ui->alpha2->setEnabled(false);
        ui->alpha3->setEnabled(false);

        // deactivate reference point

        //ui->zVectorGroup->setEnabled(false);

        //ui->z01->setEnabled(false);
        //ui->z02->setEnabled(false);
        //ui->z03->setEnabled(false);

        // deactivate reference direction (n-vector)

        ui->referencePointGroup->setEnabled(false);

        //ui->zVector1->setEnabled(false);
        //ui->zVector2->setEnabled(false);
        //ui->zVector3->setEnabled(false);

        //ui->alphaParameterGroup->hide();
    }
}

void InflowParameterWidget::on_selectExponential_clicked(bool checked)
{
    // adjust UI for exponential distribution
    if ( checked ) {

        // activate alpha parameters

        ui->alpha1->setEnabled(true);
        ui->alpha2->setEnabled(true);
        ui->alpha3->setEnabled(true);

        // activate reference point

        //ui->zVectorGroup->setEnabled(true);

        //ui->z01->setEnabled(true);
        //ui->z02->setEnabled(true);
        //ui->z03->setEnabled(true);

        // activate reference direction (n-vector)

        ui->referencePointGroup->setEnabled(true);

        //ui->zVector1->setEnabled(true);
        //ui->zVector2->setEnabled(true);
        //ui->zVector3->setEnabled(true);

        //ui->alphaParameterGroup->show();
    }
}

void InflowParameterWidget::on_resetButton_clicked()
{
    // set UI to default parameter values
    setDefaultParameters();
}

void InflowParameterWidget::on_comboBox_currentIndexChanged(int index)
{
    // this is where we get a mode
}
