#-------------------------------------------------
#
# Project created by QtCreator 2019-05-20T09:50:35
#
#-------------------------------------------------

QT       += core gui network

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = TurbulentInflowTool
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

CONFIG += c++11

NEEDED_REPO=SimCenterCommon
NEEDED_PACKAGE=Common
NEEDED_PRI=$$PWD/../$$NEEDED_REPO//$$NEEDED_PACKAGE/$$"$$NEEDED_PACKAGE".pri
!exists( $$NEEDED_PRI ) {
    message("Needed Git repo $$NEEDED_REPO not found. This project requires $$NEEDED_REPO from https://github.com/NHERI-SimCenter.")
}

include($$NEEDED_PRI)

SOURCES += \
        customizeditemmodel.cpp \
        main.cpp \
        mainwindow.cpp \
        utilWindows/dialogabout.cpp \
        utilWindows/helpwindow.cpp \
        widgets/exportwidget.cpp \
        widgets/filewidget.cpp \
        widgets/inflowparameterwidget.cpp

HEADERS += \
        customizeditemmodel.h \
        mainwindow.h \
        utilWindows/dialogabout.h \
        utilWindows/helpwindow.h \
        widgets/exportwidget.h \
        widgets/filewidget.h \
        widgets/inflowparameterwidget.h

FORMS += \
        mainwindow.ui \
        utilWindows/dialogabout.ui \
        utilWindows/helpwindow.ui \
        widgets/exportwidget.ui \
        widgets/filewidget.ui \
        widgets/inflowparameterwidget.ui

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

RESOURCES += \
    inflowtool.qrc
