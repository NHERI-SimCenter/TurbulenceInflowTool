#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QModelIndex>
#include <QItemSelection>
#include <QMap>
#include <QString>
#include <QByteArray>

class QDir;

class QStandardItemModel;

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_action_Quit_triggered();
    void on_action_About_triggered();
    void on_action_New_triggered();
    void on_action_Open_triggered();
    void on_action_Save_triggered();
    void on_actionSave_As_triggered();

    void on_btn_selectSource_clicked();

    void selectionChangedSlot(const QItemSelection & /*newSelection*/, const QItemSelection &/*oldSelection*/);

    void on_action_Documentation_triggered();

public slots:
    void fetchUFileData(bool, QDir &);

private:
    Ui::MainWindow *ui;

    QModelIndex infoItemIdx;
    QStandardItemModel *standardModel;

    QByteArray head;
    QByteArray tail;
    QMap<QString, QMap<QString, QString> * > data;
};

#endif // MAINWINDOW_H
