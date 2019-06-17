#ifndef FILEWIDGET_H
#define FILEWIDGET_H

#include <QFrame>

namespace Ui {
class FileWidget;
}

class FileWidget : public QFrame
{
    Q_OBJECT

public:
    explicit FileWidget(QWidget *parent = nullptr);
    ~FileWidget();

private slots:
    void on_sourceLocateBtn_clicked();

private:
    Ui::FileWidget *ui;

    bool validSourcePresent = false;
};

#endif // FILEWIDGET_H
