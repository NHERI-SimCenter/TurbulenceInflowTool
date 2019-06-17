#ifndef EXPORTWIDGET_H
#define EXPORTWIDGET_H

#include <QFrame>
#include <QDir>

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
    void sendParametersMap(void);

public slots:
    void setLocationAvailable(bool, QDir &);

private:
    Ui::ExportWidget *ui;

    void exportInflowParameterFile(QString);

    bool hasLocation = false;
    QDir oldLocation = QDir(".");
    QDir newLocation = QDir(".");
    QMap<QString, double> *theParameters;
};

#endif // EXPORTWIDGET_H
