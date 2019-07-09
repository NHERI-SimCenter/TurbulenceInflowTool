#ifndef INFLOWPARAMETERWIDGET_H
#define INFLOWPARAMETERWIDGET_H

#include <QFrame>
#include <QMap>

namespace Ui {
class InflowParameterWidget;
}

class InflowParameterWidget : public QFrame
{
    Q_OBJECT

public:
    explicit InflowParameterWidget(QWidget *parent = nullptr);
    ~InflowParameterWidget();
    void selectSourceLocation(void);

signals:
    void parametersReady(QMap<QString, double> &);

private slots:
    void on_RB_digitalFilter_clicked();
    void on_RB_syntheticEddie_clicked();
    void on_PHI21_valueChanged(double arg1);
    void on_PHI31_valueChanged(double arg1);
    void on_PHI32_valueChanged(double arg1);
    void setDefaultParameters();
    void on_resetButton_clicked();
    void on_modelSelectionCBX_currentIndexChanged(int index);
    void sendParameterMap(void);

private:
    void setUniformTurbulent(void);
    void setExponentialTurbulent(void);
    void refreshParameterMap(void);
    void refreshDisplay(void);

    Ui::InflowParameterWidget *ui;

    QMap<QString, double> theParameters;
    bool hasParameters = false;
};

#endif // INFLOWPARAMETERWIDGET_H
