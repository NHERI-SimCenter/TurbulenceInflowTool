#ifndef FILEWIDGET_H
#define FILEWIDGET_H

#include <QFrame>
#include <QDir>
#include <QMap>

class QStandardItemModel;

namespace Ui {
class FileWidget;
}

class FileWidget : public QFrame
{
    Q_OBJECT

public:
    explicit FileWidget(QWidget *parent = nullptr);
    ~FileWidget();
    bool fetchUFileData(QByteArray &head, QByteArray &tail, QMap<QString, QMap<QString, QString> * > &data);

signals:
    void hasValidSource(bool, QDir &);
    void boundarySelection(int);
    void sendModel(QStandardItemModel *);

public slots:
    void setBoundarySelection(int);

private slots:
    void on_sourceLocateBtn_clicked();
    void on_boundarySelection_currentIndexChanged(int index);

private:
    bool readUfile(QString);
    QStringList getLine(void);
    QMap<QString, QString> *readParameters(void);

    Ui::FileWidget *ui;

    bool validSourcePresent = false;
    QString UFilePath;
    QByteArray UFileContents = "";
    QByteArray UFileTail = "";

    QFile UFile;
    QList<QByteArray> UFileList;
    QListIterator<QByteArray> *UIter;
    QMap<QString, QMap<QString, QString> * > boundaries;
};

#endif // FILEWIDGET_H
