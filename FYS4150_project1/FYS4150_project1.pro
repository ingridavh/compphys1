TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    specialized.cpp \
    gaussian.cpp

DISTFILES += \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N10.txt \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N100.txt \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N1000.txt \
    FYS4150_project1.py

