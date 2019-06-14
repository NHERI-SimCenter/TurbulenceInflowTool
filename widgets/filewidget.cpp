#include "filewidget.h"
#include "ui_filewidget.h"

FileWidget::FileWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::FileWidget)
{
    ui->setupUi(this);
    ui->sourceSelectionBrowser->setStyleSheet("background: #eeeeee;");
}

FileWidget::~FileWidget()
{
    delete ui;
}
