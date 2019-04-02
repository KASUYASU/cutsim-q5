#-------------------------------------------------
#
# Project created by QtCreator 2018-11-16T10:45:58
#
#-------------------------------------------------

QT       += core gui xml opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

QT += concurrent

TARGET = cutsim
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

SOURCES += \
    src/app/cutsim_window.cpp \
    src/app/levelmeter.cpp \  
    src/app/lex_analyzer.cpp \
    src/app/main.cpp \
    src/app/text_area.cpp \
    src/cutsim/bbox.cpp \
    src/cutsim/cutsim.cpp \
    src/cutsim/gldata.cpp \
    src/cutsim/glwidget.cpp \
    src/cutsim/machine.cpp \
    src/cutsim/marching_cubes.cpp \
    src/cutsim/octnode.cpp \
    src/cutsim/octree.cpp \
    src/cutsim/stl.cpp \
    src/cutsim/volume.cpp \
    src/g2m/canonLine.cpp \
    src/g2m/canonMotion.cpp \
    src/g2m/canonMotionless.cpp \
    src/g2m/g2m.cpp \
    src/g2m/helicalMotion.cpp \
    src/g2m/linearMotion.cpp \
    src/g2m/machineStatus.cpp \
    src/g2m/nanotimer.cpp

HEADERS += \
    src/app/cutsim_app.hpp \
    src/app/cutsim_def.hpp \
    src/app/cutsim_window.hpp \
    src/app/levelmeter.hpp \
    src/app/lex_analyzer.hpp \
    src/app/text_area.hpp \
    src/cutsim/bbox.hpp \
    src/cutsim/cube_wireframe.hpp \
    src/cutsim/cutsim.hpp \
    src/cutsim/facet.hpp \
    src/cutsim/gldata.hpp \
    src/cutsim/glwidget.hpp \
    src/cutsim/glvertex.hpp \
    src/cutsim/isosurface.hpp \
    src/cutsim/machine.hpp \
    src/cutsim/marching_cubes.hpp \
    src/cutsim/octnode.hpp \
    src/cutsim/octree.hpp \
    src/cutsim/stl.hpp \
    src/cutsim/volume.hpp \
    src/g2m/canonLine.hpp \
    src/g2m/canonMotion.hpp \
    src/g2m/canonMotionless.hpp \
    src/g2m/g2m.hpp \
    src/g2m/gplayer.hpp \
    src/g2m/helicalMotion.hpp \
    src/g2m/linearMotion.hpp \
    src/g2m/machineStatus.hpp \
    src/g2m/nanotimer.hpp \
    src/g2m/point.hpp \
    src/app/version_string.hpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

QMAKE_CFLAGS += -march=native -ftree-vectorize -O3 -g
QMAKE_CXXFLAGS += -march=native -ftree-vectorize -O3 -g
#QMAKE_CFLAGS += -march=native -ftree-vectorize -O0 -g
#QMAKE_CXXFLAGS += -march=native -ftree-vectorize -O0 -g

LIBS *= -L$${LIB_DIR} -lQGLViewer -lGLU -lgomp
LIB_NAME = QGLViewer-qt5
LIBS *= -L$${LIB_DIR} -l$${LIB_NAME}
