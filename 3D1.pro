TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
    TP_OPENGL.cpp \
    CVector.cpp \
    CPoint.cpp \
    glad.c \
    glew.c

DISTFILES += \
    TDTP_Intro3DOpenGL.pdf \
    main2 \
    3D1.pro.user

HEADERS += \
    CPoint.h \
    CVector.h

LIBS += -lglut
LIBS += -lGL
LIBS += -lGLU
LIBS += -ldl
