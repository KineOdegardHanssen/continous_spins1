TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp\
    bond.cpp \
    site.cpp \
    lattice.cpp \
    montecarlo.cpp \
    gaussiandeviate.cpp

HEADERS += \
    bond.h \
    site.h \
    lattice.h \
    bond.h \
    lattice.h \
    site.h \
    montecarlo.h \
    gaussiandeviate.h

LIBS+= -lfftw3

QMAKE_CXXFLAGS += -Wall -std=c++0x
INCLUDEPATH += "/home/ubu/Downloads/fftw-3.3.6-pl1"
INCLUDEPATH += "/usr/share/doc/libfftw3-3"
