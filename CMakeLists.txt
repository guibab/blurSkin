cmake_minimum_required(VERSION 2.6)
project(brSkinBrush)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
SET(MAYA_VERSION 2020 CACHE STRING "Maya version number")

SET(SOURCE_FILES 
    "src/blurSkinPlugin.cpp" 
    "src/blurSkinCmd.cpp" 
    "src/blurSkinCmd.h"
    "src/blurSkinEdit.cpp"
    "src/blurSkinEdit.h"
    "src/functions.cpp"    
    "src/functions.h"
    "src/pointsDisplay.cpp"
    "src/pointsDisplay.h"
)

FIND_PACKAGE(OpenGL REQUIRED)
FIND_PACKAGE(Maya REQUIRED)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../output)
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT ${PROJECT_NAME})

include_directories(${MAYA_INCLUDE_DIR} maya lib ${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS})
link_directories(${MAYA_LIBRARY_DIR})
add_library(${PROJECT_NAME} SHARED ${SOURCE_FILES})
target_link_libraries(${PROJECT_NAME} ${MAYA_LIBRARIES} ${OPENGL_LIBRARIES})

MAYA_PLUGIN(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} ${MAYA_TARGET_TYPE} DESTINATION plug-ins/win64/${MAYA_VERSION})
#install(CODE "execute_process(COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/scripts/buildModFile.py)")
