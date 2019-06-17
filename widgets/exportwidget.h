#ifndef EXPORTWIDGET_H
#define EXPORTWIDGET_H

#include <QFrame>
#include <QDir>
#include <QMap>

namespace Ui {
class ExportWidget;
}

class ExportWidget : public QFrame
{
    Q_OBJECT

public:
    explicit ExportWidget(QWidget *parent = nullptr);
    ~ExportWidget();

signals:
    void sendParameterMap(void);

public slots:
    void setLocationAvailable(bool, QDir &);
    void setParameterMap(QMap<QString, double> &);

private slots:
    void on_btn_export_clicked();

private:
    Ui::ExportWidget *ui;

    void exportInflowParameterFile(QString);

    bool hasLocation = false;
    bool hasParameters = false;
    QDir oldLocation = QDir(".");
    QDir newLocation = QDir(".");
    QMap<QString, double> theParameters;
};

#endif // EXPORTWIDGET_H
