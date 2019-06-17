#include "helpwindow.h"
#include "ui_helpwindow.h"

#include <QGuiApplication>
#include <QSize>
#include <QScreen>

#include <QUrl>

HelpWindow::HelpWindow(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::HelpWindow)
{
    ui->setupUi(this);

    //
    // adjust size of application window to the available display
    //
    QSize rec = QGuiApplication::primaryScreen()->size();
    int height = this->height()<0.5*rec.height()?0.5*rec.height():this->height();
    int width  = this->width()<0.5*rec.width()?0.5*rec.width():this->width();
    this->resize(width, height);
}

HelpWindow::~HelpWindow()
{
    delete ui;
}

void HelpWindow::on_pushButton_clicked()
{
    this->close();
}
