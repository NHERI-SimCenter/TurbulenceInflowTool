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

private:
    Ui::FileWidget *ui;
};

#endif // FILEWIDGET_H
