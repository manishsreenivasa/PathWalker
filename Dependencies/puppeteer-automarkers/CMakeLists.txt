CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

# Needed for UnitTest++
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake )

# Perform the proper linking
SET (CMAKE_SKIP_BUILD_RPATH FALSE)
SET (CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
SET (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

SET( QT_USE_QTXML TRUE )
SET( QT_USE_QTOPENGL TRUE )
SET( QT_WRAP_CPP TRUE )
FIND_PACKAGE (Qt4 REQUIRED QtCore QtGui)
INCLUDE (${QT_USE_FILE})

FIND_PACKAGE (Lua REQUIRED)
FIND_PACKAGE (OpenGL)
FIND_PACKAGE (RBDL COMPONENTS LuaModel REQUIRED)
FIND_PACKAGE (Eigen3 REQUIRED)

ADD_DEFINITIONS (${QT_DEFINITIONS})

QT4_WRAP_UI (PuppeteerMainWindow_UIS_H ${CMAKE_CURRENT_SOURCE_DIR}/ui/PuppeteerMainWindow.ui)

MESSAGE (STATUS "out ${PuppeteerMainWindow_UIS_H}")

SET ( SCENEGL_SRCS
	src/Camera.cc
	src/Scene.cc
	src/MeshVBO.cc
	src/Shader.cc
	src/Model.cc
	src/MarkerData.cc
	src/Animation.cc
	src/ModelFitter.cc
	src/Scripting.cc
	)

SET ( QtGLBaseApp_SRCS
  src/main.cc
  src/GLWidget.cc
  src/PuppeteerApp.cc
  glew/src/glew.c
  ## CHART STUFF
  src/vtkChart/chartXY.cc
  )

QT4_WRAP_CPP ( QtGLBaseApp_MOC_SRCS
  src/PuppeteerApp.h
  src/GLWidget.h
  )

ADD_LIBRARY ( SceneGL STATIC
  ${SCENEGL_SRCS}
	${PuppeteerMainWindow_UIS_H} 
	)

ADD_LIBRARY ( glew STATIC
  glew/src/glew.c
  )

ADD_EXECUTABLE ( puppeteer
  ${QtGLBaseApp_SRCS}
  ${QtGLBaseApp_MOC_SRCS}
	${PuppeteerMainWindow_UIS_H}
  )

ADD_EXECUTABLE ( fit_motion_automarkers
  src/fit_motion_automarkers.cc
  )

INCLUDE_DIRECTORIES (
  ${QT_INCLUDE_DIR}
  glew/include
  )

SET ( BUILD_SOURCE_DIRECTORY ${CMAKE_SOURCE_DIR} )
SET ( BUILD_INSTALL_DIRECTORY ${CMAKE_INSTALL_PREFIX} )
CONFIGURE_FILE( ${CMAKE_SOURCE_DIR}/src/config.h.cmake ${CMAKE_BINARY_DIR}/src/config.h )

INCLUDE_DIRECTORIES (
	${QT_INCLUDE_DIR}
	glew/include
	${QT_INCLUDE_DIR}/QtOpenGL
	${CMAKE_CURRENT_BINARY_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}/src
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/QtPropertyBrowser/src
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/c3dfile/src
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/luatables
	${CMAKE_CURRENT_SOURCE_DIR}/vendor/tclap/include
	${CMAKE_BINARY_DIR}/src
	${RBDL_INCLUDE_DIR}
	${EIGEN3_INCLUDE_DIR}
	${LUA_INCLUDE_DIR}
	)

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
  LINKER_LANGUAGE CXX
  INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib" 
  )

SUBDIRS (
	tests
	vendor/QtPropertyBrowser
	vendor/c3dfile
	vendor/luatables
	)

TARGET_LINK_LIBRARIES ( SceneGL
	glew
	luatables
	${LUA_LIBRARY}
	${OPENGL_LIBRARIES}
	${RBDL_LIBRARY}
	${RBDL_LuaModel_LIBRARY}
	c3dfile
	)

if(VTK_LIBRARIES)
  if(${VTK_VERSION} VERSION_LESS "6")
    target_link_libraries(puppeteer 
      SceneGL
      ${QT_LIBRARIES}
      QtPropertyBrowser
      ${OPENGL_LIBRARIES}
      ${Boost_SYSTEM_LIBRARY}
			${VTK_LIBRARIES} 
			QVTK
      )
  else()
    target_link_libraries(puppeteer 
      SceneGL
      ${QT_LIBRARIES}
      QtPropertyBrowser
      ${Boost_SYSTEM_LIBRARY}
			${VTK_LIBRARIES}
      )
  endif()
else()
  target_link_libraries(
    puppeteer
    SceneGL
    ${QT_LIBRARIES}
    QtPropertyBrowser
    ${Boost_SYSTEM_LIBRARY}
    vtkHybrid 
    QVTK 
    vtkViews 
		${QT_LIBRARIES}
    )
endif()

TARGET_LINK_LIBRARIES ( fit_motion_automarkers
  SceneGL
	-lpthread # fix_nvidia_linking_new_ubuntu
  )

# Installation
INSTALL (TARGETS puppeteer fit_motion_automarkers
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  )

# Packaging
SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_SECTION "science")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Martin Felis <martin.felis@iwr.uni-heidelberg.de>")
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libqt4-opengl (>= 4:4.8), libglu1-mesa
(>=8.0), rbdl (> 2:2.2)")
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Puppeteer - Motion Capture Mapping Tool")
SET(CPACK_PACKAGE_VENDOR "Martin Felis")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
SET(CPACK_PACKAGE_VERSION_MAJOR ${MESHUP_VERSION_MAJOR})
SET(CPACK_PACKAGE_VERSION_MINOR ${MESHUP_VERSION_MINOR})
SET(CPACK_PACKAGE_VERSION_PATCH ${MESHUP_VERSION_PATCH})
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CPACK_PACKAGE ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}")
SET(CPACK_PACKAGE_FILE_NAME "MeshUp-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
IF(UNIX)
  SET(CPACK_STRIP_FILES "bin/puppeteer")
  SET(CPACK_SOURCE_STRIP_FILES "")
ENDIF(UNIX)
SET(CPACK_PACKAGE_EXECUTABLES "puppeteer;Puppeteer")

INCLUDE(CPack)
