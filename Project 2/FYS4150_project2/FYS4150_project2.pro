TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -larmadillo

SOURCES += main.cpp \
    jacobi_rotation.cpp \
    schrodinger.cpp

HEADERS += \
    jacobi_rotation.h \
    schrodinger.h

DISTFILES += \
    Benchmarks.txt

