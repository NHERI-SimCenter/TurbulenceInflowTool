#ifndef EXPORTWIDGET_H
#define EXPORTWIDGET_H

#include <QFrame>
#include <QDir>
#include <QMap>

class QStandardItemModel;

namespace Ui {
class ExportWidget;
}

class ExportWidget : public QFrame
{
    Q_OBJECT

public:
    explicit ExportWidget(QWidget *parent = nullptr);
    ~ExportWidget();
    void setUFileData(QByteArray &head, QByteArray &tail, QMap<QString, QMap<QString, QString> * > &data);

signals:
    void sendParameterMap(void);
    void boundarySelection(int);

public slots:
    void setLocationAvailable(bool, QDir &);
    void setParameterMap(QMap<QString, double> &);
    void setBoundarySelection(int);
    void setModel(QStandardItemModel *);

private slots:
    void on_btn_export_clicked();
    void on_boundarySelection_currentIndexChanged(int index);

private:
    Ui::ExportWidget *ui;

    void exportInflowParameterFile(QString);
    void exportUFile(QString);
    void clearBoundaryMap(void);

    bool hasLocation = false;
    bool hasParameters = false;
    QDir oldLocation = QDir(".");
    QDir newLocation = QDir(".");
    QMap<QString, double> theParameters;

    bool validSourcePresent = false;
    QString UFilePath;

    QFile UFile;
    QList<QByteArray> UFileList;
    QListIterator<QByteArray> *UIter;

    QByteArray UFileHead = "";
    QByteArray UFileTail = "";
    QMap<QString, QMap<QString, QString> * > boundaries;
};

#endif // EXPORTWIDGET_H
