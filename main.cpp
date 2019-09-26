#include "mainwindow.h"
#include <QApplication>
#include <QFile>
#include <QGuiApplication>

QString openStyleFiles()
{
    QString ret;
    QFile mainStyleFile(":/resources/styleSheets/TInFStyles.qss");

#ifdef Q_OS_WIN
    QFile appendedStyle(":/resources/styleSheets/TInFWin.qss");
#endif

#ifdef Q_OS_MACOS
    QFile appendedStyle(":/resources/styleSheets/TInFMac.qss");
#endif

#ifdef Q_OS_LINUX
    QFile appendedStyle(":/resources/styleSheets/TInFLinux.qss");
#endif

    if (!mainStyleFile.open(QFile::ReadOnly))
    {
        return ret;
    }

    if (!appendedStyle.open(QFile::ReadOnly))
    {
        return ret;
    }

    ret = ret.append(mainStyleFile.readAll());
    ret = ret.append(appendedStyle.readAll());

    mainStyleFile.close();
    appendedStyle.close();

    return ret;
}

int main(int argc, char *argv[])
{
    QGuiApplication::setAttribute(Qt::AA_EnableHighDpiScaling);

    QApplication app(argc, argv);
    MainWindow w;
    QApplication::setWindowIcon(QIcon(":/resources/NHERI-TInF-icon.icns"));
    w.show();

    app.setStyleSheet(openStyleFiles());

    return app.exec();
}
