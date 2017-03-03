CMAKE_MINIMUM_REQUIRED (VERSION 2.6)

CMAKE_POLICY(SET CMP0048 NEW)
CMAKE_POLICY(SET CMP0040 NEW)
 
SET ( RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_MAJOR 1 )
SET ( RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_MINOR 0 )
SET ( RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_PATCH 0 )

SET ( RBDL_ADDON_CUSTOMFORCES_VERSION 
	${RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_MAJOR}.${RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_MINOR}.${RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION_PATCH}
)


PROJECT (RBDL_CUSTOMFORCES_TESTS VERSION ${RBDL_ADDON_CUSTOMFORCES_VERSION})
#SET (PROJECT_VERSION ${RBDL_ADDON_CUSTOMFORCES_TESTS_VERSION})

# Needed for UnitTest++
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../CMake )

# Look for unittest++
FIND_PACKAGE (UnitTest++ REQUIRED)
INCLUDE_DIRECTORIES (${UNITTEST++_INCLUDE_DIR})

SET ( CUSTOMFORCES_TESTS_SRCS
	testAnderson2007TorqueMuscle.cc
	testNeuroSwingTorsionSpring.cc
    ../customforces.h
    ../Anderson2007TorqueMuscle.h
    ../Anderson2007TorqueMuscle.cc
    ../NeuroSwingTorsionSpring.h
    ../NeuroSwingTorsionSpring.cc
    ../csvtools.h
    ../csvtools.cc
    ../../geometry/geometry.h
    ../../geometry/SegmentedQuinticBezierToolkit.h
    ../../geometry/SmoothSegmentedFunction.h
    ../../geometry/SmoothSegmentedFunctionFactory.h
    ../../geometry/SegmentedQuinticBezierToolkit.cc
    ../../geometry/SmoothSegmentedFunction.cc
    ../../geometry/SmoothSegmentedFunctionFactory.cc
	)

INCLUDE_DIRECTORIES ( ../ ../../geometry )

SET_TARGET_PROPERTIES ( ${PROJECT_EXECUTABLES} PROPERTIES
  LINKER_LANGUAGE CXX
)

ADD_EXECUTABLE ( rbdl_customforces_tests ${CUSTOMFORCES_TESTS_SRCS} )

SET_TARGET_PROPERTIES ( rbdl_customforces_tests PROPERTIES
	LINKER_LANGUAGE CXX
	OUTPUT_NAME runCustomForcesTests
	)

SET (RBDL_LIBRARY rbdl)
IF (RBDL_BUILD_STATIC)
	SET (RBDL_LIBRARY rbdl-static)
ENDIF (RBDL_BUILD_STATIC)

TARGET_LINK_LIBRARIES ( rbdl_customforces_tests
		${UNITTEST++_LIBRARY}
		${RBDL_LIBRARY}
	)
	
OPTION (RUN_AUTOMATIC_TESTS "Perform automatic tests after compilation?" OFF)

IF (RUN_AUTOMATIC_TESTS)
ADD_CUSTOM_COMMAND (TARGET rbdl_customforces_tests
	POST_BUILD
	COMMAND ./runCustomForcesTests
	COMMENT "Running automated addon customforces tests..."
	)
ENDIF (RUN_AUTOMATIC_TESTS)
