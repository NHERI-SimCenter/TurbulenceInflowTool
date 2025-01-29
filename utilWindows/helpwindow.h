#ifndef HELPWINDOW_H
#define HELPWINDOW_H

#include <QFrame>

namespace Ui {
class HelpWindow;
}

class HelpWindow : public QFrame
{
    Q_OBJECT

public:
    explicit HelpWindow(QWidget *parent = nullptr);
    ~HelpWindow();

private slots:
    void on_pushButton_clicked();

private:
    Ui::HelpWindow *ui;
};

#endif // HELPWINDOW_H
