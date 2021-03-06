TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    specialized.cpp \
    gaussian.cpp \
    my_lu.cpp

DISTFILES += \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N10.txt \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N100.txt \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_N1000.txt \
    FYS4150_project1.py \
    ../build-FYS4150_project1-Desktop-Debug/p1_result_s_N10.txt

HEADERS += \
    gaussian.h \
    specialized.h \
    my_lu.h


unix:!macx: LIBS += -larmadillo
