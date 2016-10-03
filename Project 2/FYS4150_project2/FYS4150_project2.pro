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
    Benchmarks.txt \
    ../build-FYS4150_project2-Desktop-Debug/eigenvec_omega=0.01.txt \
    ../build-FYS4150_project2-Desktop-Debug/eigenvec_omega=0.5.txt \
    ../build-FYS4150_project2-Desktop-Debug/eigenvec_omega=1.txt \
    ../build-FYS4150_project2-Desktop-Debug/eigenvec_omega=5.txt \
    ../FYS4150_project2.py

