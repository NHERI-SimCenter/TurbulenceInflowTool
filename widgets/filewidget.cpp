#include "filewidget.h"
#include "ui_filewidget.h"

#include <QString>
#include <QDir>
#include <QFileDialog>

#include <QDebug>

FileWidget::FileWidget(QWidget *parent) :
    QFrame(parent),
    ui(new Ui::FileWidget)
{
    ui->setupUi(this);
}

FileWidget::~FileWidget()
{
    delete ui;
}

void FileWidget::on_sourceLocateBtn_clicked()
{
    // select the openFOAM input source tree
    QFileDialog *dlg = new QFileDialog();
    dlg->setReadOnly(true);
    dlg->setFileMode(QFileDialog::DirectoryOnly);
    dlg->exec();
    QDir fileTreeLocation = dlg->directory();

    ui->sourceLocationDisplay->setText(fileTreeLocation.path());

    QStringList folders = fileTreeLocation.entryList(QStringList(),QDir::Dirs);
    int stack = folders.length();

    if (folders.contains("0") && folders.contains("constant") && folders.contains("system") ) {
        ui->sourceLocationDisplay->setStyleSheet("color: #000000;");
        validSourcePresent = true;
    }
    else {
        ui->sourceLocationDisplay->setStyleSheet("color: #ff0000;");
        validSourcePresent = false;
    }

    delete dlg;

    if (validSourcePresent) {
        // parse files for available boundaries

        // let main application know that source is available
        emit hasValidSource(true, fileTreeLocation);
    }
    else {
        // user not ready to proceed
        QDir thisDir(".");
        emit hasValidSource(false, thisDir);
    }
}
