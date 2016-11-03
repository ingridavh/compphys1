TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    spins.cpp \
    spinsystem.cpp

DISTFILES += \
    benchmarks \
    benchmarks.txt \
    log.txt

HEADERS += \
    spinsystem.h

