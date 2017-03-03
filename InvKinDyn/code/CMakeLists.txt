CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT (InvKinDyn CXX)

# We need to add the project source path to the CMake module path so that
# the FindRBDL.cmake script can be found.
LIST( APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} )

# Search for the RBDL include directory and library
FIND_PACKAGE (RBDL 2.5.0 REQUIRED)
FIND_PACKAGE (Eigen3 3.0.0 REQUIRED)
FIND_PACKAGE (BTK REQUIRED)

# Add the include directory to the include paths
INCLUDE(${BTK_USE_FILE})
INCLUDE_DIRECTORIES ( ${RBDL_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${BTK_INCLUDE_DIR})

# Create an executable
ADD_EXECUTABLE (compute_id compute_id.cc)

# And link the library against the executable
TARGET_LINK_LIBRARIES ( compute_id ${RBDL_LIBRARIES} ${BTK_LIBRARIES})

