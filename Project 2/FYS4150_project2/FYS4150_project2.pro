TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo

SOURCES += main.cpp \
    jacobi.cpp \
    jacobi_rotation.cpp

HEADERS += \
    jacobi.h \
    jacobi_rotation.h

DISTFILES += \
    Benchmarks.txt

