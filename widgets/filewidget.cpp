#include "filewidget.h"
#include "ui_filewidget.h"

#include <QString>
#include <QDir>
#include <QFileDialog>
#include <QList>
#include <QListIterator>
#include <QStandardItemModel>

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
    UFileHead = "";
    UFileTail = "";

    // select the openFOAM input source tree
    QFileDialog *dlg = new QFileDialog();
    dlg->setReadOnly(true);
    dlg->setFileMode(QFileDialog::DirectoryOnly);
    dlg->exec();
    QDir fileTreeLocation = dlg->directory();

    ui->sourceLocationDisplay->setText(fileTreeLocation.path());

    // if selection made, open it and check valid dir .. if valid
    // set fileTreeLocation to point to it
    QStringList filesSelected;
    filesSelected = dlg->selectedFiles();
    if (filesSelected.size() == 1) {
        QString dirName = filesSelected.at(0);
        QDir dir(dirName);
        if (dir.exists())
            fileTreeLocation = dir;
    }

    QStringList folders = fileTreeLocation.entryList(QStringList(),QDir::Dirs);

    //QStringList folders = fileTreeLocation.entryList(filesSelected);
    int stack = folders.length();

    if (folders.contains("0") && folders.contains("constant")  && folders.contains("system")) {
        //
        // look for U file
        //

        QDir UDir = fileTreeLocation;
        UDir.cd("0");
        UFilePath = UDir.filePath("U");

        //
        // look for controlDict file
        //

        QDir CDictDir = fileTreeLocation;
        CDictDir.cd("system");
        ControlDictPath = CDictDir.filePath("controlDict");

        validSourcePresent = readUfile(UFilePath) && readControlDict(ControlDictPath);

        if (validSourcePresent)
        { ui->sourceLocationDisplay->setStyleSheet("color: #000000;"); }
        else
        { ui->sourceLocationDisplay->setStyleSheet("color: #FFA500;"); }
    }
    else {
        //
        // this is not a valid OpenFOAM folder
        //
        UFileContents = "";

        ui->sourceLocationDisplay->setStyleSheet("color: #ff0000;");
        validSourcePresent = false;
    }

    delete dlg;

    if (validSourcePresent) {
        // parse files for available boundaries
        QStringList boundaryList;

        UFileList = UFileContents.split('\n');
        UIter = new QListIterator<QByteArray>(UFileList);

        // read till boundaryField keyword
        while (UIter->hasNext())
        {
            QByteArray line = UIter->next();
            UFileHead.append(line);
            UFileHead.append('\n');
            if (line.contains("boundaryField")) {
                while ( (!line.contains('{')) && UIter->hasNext()) {
                    line = UIter->next();
                    UFileHead.append(line);
                    UFileHead.append('\n');
                }
                break;
            }
        }

        // parse for boundary patches
        while (UIter->hasNext())
        {
            QStringList list;

            if (this->getLine(list))
            {
                // skip empty lines
                if (list.length() == 0) continue;

                // terminate if done with boundaryFields section
                if (list[0] == '}') {
                    UFileTail.append("}\n");
                    break;
                }

                // read and store the boundary item
                boundaryList.append(list[0]);
                boundaries.insert(list[0], this->readParameters());
            }
        }

        // collect the remainder of the file
        while (UIter->hasNext())
        {
            QByteArray line = UIter->next();
            UFileTail.append(line);
            UFileTail.append('\n');
        }

        // let main application know that source is available
        emit hasValidSource(true, fileTreeLocation);

        QStandardItemModel *theModel= new QStandardItemModel();
        foreach(QString s, boundaryList)
        {
            theModel->appendRow(new QStandardItem(s));
        }
        ui->boundarySelection->setModel(theModel);
        emit sendModel(theModel);
    }
    else {
        // user not ready to proceed
        QDir thisDir(".");
        emit hasValidSource(false, thisDir);
    }

#if 0
    foreach (QString k, boundaries.keys())
    {
        qDebug() << k << ": " << *(boundaries.value(k));
    }
#endif
}

bool FileWidget::readUfile(QString filename)
{
    QFile UFile(filename);

    if (UFile.exists()) {
        //
        // U file exists
        //
        UFile.open(QFile::ReadOnly);
        UFileContents = UFile.readAll();
        UFile.close();

        return true;
    }
    else {
        //
        // U file missing
        //
        UFileContents = "";

        return false;
    }
}

bool FileWidget::readControlDict(QString filename)
{
    QFile CDictFile(filename);

    if (CDictFile.exists()) {
        //
        // controlDict file exists
        //
        CDictFile.open(QFile::ReadOnly);
        CDictContents = CDictFile.readAll();
        CDictFile.close();

        return true;
    }
    else {
        //
        // controlDict file missing
        //
        CDictContents = "";

        return false;
    }
}

bool FileWidget::getLine(QStringList &reply)
{
    bool hasLine = false;
    QByteArray lineString = "";

    while (UIter->hasNext() && (!hasLine))
    {
        QByteArray line = UIter->next().simplified();
        if (qstrncmp(line,"//",2) == 0) continue;
        if (qstrncmp(line, "#",1) == 0) continue;
        if (line.contains('{')) {
            hasLine = true;
            break;
        }
        lineString += line;
        if (line.contains('}')) {
            hasLine = true;
            break;
        }
        if (line.contains(';')) {
            int idx = lineString.indexOf(';');
            lineString.truncate(idx);
            hasLine = true;
            break;
        }
    }

    reply.clear();

    if (hasLine)
    {
        QByteArrayList reply0 = lineString.simplified().split(' ');

        foreach (QByteArray item, reply0)
        {
            reply.append(item);
        }
    }

    return hasLine;
}

QMap<QString, QString> *FileWidget::readParameters(void)
{
    QMap<QString, QString> *params = new QMap<QString, QString>();

    QStringList items;

    while ( this->getLine(items) ) {
        if (items[0] == '}') break;

        if (items.length() > 0 ) {
            QString key = items[0];
            items.removeFirst();
            QString value = items.join(" ");
            params->insert(key, value);
        }
    }

    return params;
}

bool FileWidget::fetchUFileData(QByteArray &head, QByteArray &tail, QMap<QString, QMap<QString, QString> * > &data)
{
    head = UFileHead;
    tail = UFileTail;
    data = boundaries;

    return (UFileContents.length()>0);
}

void FileWidget::setBoundarySelection(int index)
{
    ui->boundarySelection->setCurrentIndex(index);
}

void FileWidget::on_boundarySelection_currentIndexChanged(int index)
{
    emit boundarySelection(index);
}
