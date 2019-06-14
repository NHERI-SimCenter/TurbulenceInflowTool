#ifndef EXPORTWIDGET_H
#define EXPORTWIDGET_H

#include <QFrame>

namespace Ui {
class ExportWidget;
}

class ExportWidget : public QFrame
{
    Q_OBJECT

public:
    explicit ExportWidget(QWidget *parent = nullptr);
    ~ExportWidget();

private:
    Ui::ExportWidget *ui;
};

#endif // EXPORTWIDGET_H
