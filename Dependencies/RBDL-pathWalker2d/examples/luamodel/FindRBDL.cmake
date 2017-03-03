# Searches for RBDL includes and library files
#
# Sets the variables
#   RBDL_FOUND
#   RBDL_INCLUDE_DIR
#   RBDL_LIBRARY

SET (RBDL_FOUND FALSE)
SET (RBDL_LuaModel_FOUND FALSE)
SET (RBDL_URDFReader_FOUND FALSE)
SET (RBDL_GEOMETRY_FOUND FALSE)
SET (RBDL_CUSTOMFORCES_FOUND FALSE)

# If you have multiple versions of RBDL installed on your
# machine, or if RBDL is not being found (because it isn't 
# in the usual search directions) you can force CMake to
# use a specific directory:
#
# Simply set the CUSTOM_RBDL_PATH to the correct folder. 
# For example:
#
#  SET(CUSTOM_RBDL_PATH home/user/dev/rbdlCustomFolderInstall )
#
# If you do not set CUSTOM_RBDL_PATH to a directory then
# CMake will look in a series of folders for RBDL
#    $ENV{HOME}/local
#    $ENV{RBDL_PATH}
#    /usr/local
#    /usr
#
SET(CUSTOM_RBDL_PATH /home/mjhmilla/dev/rbdlAddonGeometryCustomForces-install)



IF(CUSTOM_RBDL_PATH)
  FIND_PATH (RBDL_INCLUDE_DIR rbdl/rbdl.h
    HINTS
    ${CUSTOM_RBDL_PATH}/include
    NO_DEFAULT_PATH
    )

  FIND_LIBRARY (RBDL_LIBRARY NAMES rbdl
    PATHS
    ${CUSTOM_RBDL_PATH}/lib
    NO_DEFAULT_PATH
    )

  FIND_PATH (RBDL_LuaModel_INCLUDE_DIR rbdl/addons/luamodel/luamodel.h
    HINTS
    ${CUSTOM_RBDL_PATH}/include
    NO_DEFAULT_PATH
    )

  FIND_LIBRARY (RBDL_LuaModel_LIBRARY NAMES rbdl_luamodel
    PATHS
    ${CUSTOM_RBDL_PATH}/lib
    NO_DEFAULT_PATH
    )

  FIND_PATH (RBDL_URDFReader_INCLUDE_DIR rbdl/addons/urdfreader/urdfreader.h
    HINTS
    ${CUSTOM_RBDL_PATH}/include
    NO_DEFAULT_PATH  
    )

  FIND_LIBRARY (RBDL_URDFReader_LIBRARY NAMES rbdl_urdfreader
    PATHS
    ${CUSTOM_RBDL_PATH}/lib
    NO_DEFAULT_PATH  
    )


  FIND_PATH (RBDL_CUSTOMFORCES_INCLUDE_DIR rbdl/addons/customforces/customforces.h
    HINTS
    ${CUSTOM_RBDL_PATH}/include
    NO_DEFAULT_PATH
    )

  FIND_LIBRARY (RBDL_CUSTOMFORCES_LIBRARY NAMES rbdl_customforces  
    PATHS
    ${CUSTOM_RBDL_PATH}/lib
    NO_DEFAULT_PATH
    )

  FIND_PATH (RBDL_GEOMETRY_INCLUDE_DIR rbdl/addons/geometry/geometry.h
    HINTS
    ${CUSTOM_RBDL_PATH}/include
    NO_DEFAULT_PATH  
    )

  FIND_LIBRARY (RBDL_GEOMETRY_LIBRARY NAMES rbdl_geometry  
    PATHS
    ${CUSTOM_RBDL_PATH}/lib
    NO_DEFAULT_PATH  
    )

ELSE(CUSTOM_RBDL_PATH)
  FIND_PATH (RBDL_INCLUDE_DIR rbdl/rbdl.h
    HINTS
    $ENV{HOME}/local/include
    $ENV{RBDL_PATH}/src
    $ENV{RBDL_PATH}/include
    $ENV{RBDL_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    )

  FIND_LIBRARY (RBDL_LIBRARY NAMES rbdl
    PATHS
    $ENV{HOME}/local/lib
    $ENV{HOME}/local/lib/x86_64-linux-gnu
    $ENV{RBDL_PATH}/lib
    $ENV{RBDL_LIBRARY_PATH}
    /usr/local/lib
    /usr/local/lib/x86_64-linux-gnu
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    )

  FIND_PATH (RBDL_LuaModel_INCLUDE_DIR rbdl/addons/luamodel/luamodel.h
    HINTS
    $ENV{HOME}/local/include
    $ENV{RBDL_PATH}/src
    $ENV{RBDL_PATH}/include
    $ENV{RBDL_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    )

  FIND_LIBRARY (RBDL_LuaModel_LIBRARY NAMES rbdl_luamodel
    PATHS
    $ENV{HOME}/local/lib
    $ENV{HOME}/local/lib/x86_64-linux-gnu
    $ENV{RBDL_PATH}
    $ENV{RBDL_LIBRARY_PATH}
    /usr/local/lib
    /usr/local/lib/x86_64-linux-gnu
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    )

  FIND_PATH (RBDL_URDFReader_INCLUDE_DIR rbdl/addons/urdfreader/urdfreader.h
    HINTS
    $ENV{HOME}/local/include
    $ENV{RBDL_PATH}/src
    $ENV{RBDL_PATH}/include
    $ENV{RBDL_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    )

  FIND_LIBRARY (RBDL_URDFReader_LIBRARY NAMES rbdl_urdfreader
    PATHS
    $ENV{HOME}/local/lib
    $ENV{HOME}/local/lib/x86_64-linux-gnu
    $ENV{RBDL_PATH}
    $ENV{RBDL_LIBRARY_PATH}
    /usr/local/lib
    /usr/local/lib/x86_64-linux-gnu
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    )


  FIND_PATH (RBDL_CUSTOMFORCES_INCLUDE_DIR rbdl/addons/customforces/customforces.h
    HINTS
    $ENV{HOME}/local/include
    $ENV{RBDL_PATH}/src
    $ENV{RBDL_PATH}/include
    $ENV{RBDL_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    )

  FIND_LIBRARY (RBDL_CUSTOMFORCES_LIBRARY NAMES rbdl_customforces  
    PATHS
    $ENV{HOME}/local/lib
    $ENV{HOME}/local/lib/x86_64-linux-gnu
    $ENV{RBDL_PATH}
    $ENV{RBDL_LIBRARY_PATH}
    /usr/local/lib
    /usr/local/lib/x86_64-linux-gnu
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    )

  FIND_PATH (RBDL_GEOMETRY_INCLUDE_DIR rbdl/addons/geometry/geometry.h
    HINTS
    $ENV{HOME}/local/include
    $ENV{RBDL_PATH}/src
    $ENV{RBDL_PATH}/include
    $ENV{RBDL_INCLUDE_PATH}
    /usr/local/include
    /usr/include
    )

  FIND_LIBRARY (RBDL_GEOMETRY_LIBRARY NAMES rbdl_geometry  
    PATHS
    $ENV{HOME}/local/lib
    $ENV{HOME}/local/lib/x86_64-linux-gnu
    $ENV{RBDL_PATH}
    $ENV{RBDL_LIBRARY_PATH}
    /usr/local/lib
    /usr/local/lib/x86_64-linux-gnu
    /usr/lib
    /usr/lib/x86_64-linux-gnu
    )  
ENDIF(CUSTOM_RBDL_PATH)

IF (RBDL_INCLUDE_DIR AND RBDL_LIBRARY)
  SET (RBDL_FOUND TRUE)
ELSE(RBDL_INCLUDE_DIR AND RBDL_LIBRARY)
  IF(RBDL_FIND_REQUIRED)
    MESSAGE (SEND_ERROR " Could not find RBDL.")
    MESSAGE (SEND_ERROR " Try setting CUSTOM_RBDL_PATH in FindRBDL.cmake force CMake to use the desired directory.")
  ELSE(RBDL_FIND_REQUIRED)
    MESSAGE (STATUS " Could not find RBDL.")
    MESSAGE (STATUS " Try setting CUSTOM_RBDL_PATH in FindRBDL.cmake force CMake to use the desired directory.")
  ENDIF(RBDL_FIND_REQUIRED)

ENDIF (RBDL_INCLUDE_DIR AND RBDL_LIBRARY)


IF (RBDL_LuaModel_LIBRARY AND RBDL_LuaModel_INCLUDE_DIR)
  SET (RBDL_LuaModel_FOUND TRUE)
ELSE (RBDL_LuaModel_LIBRARY AND RBDL_LuaModel_INCLUDE_DIR)
  MESSAGE(SEND_ERROR "Could not find RDBL ADDONS LUAMODEL")
ENDIF (RBDL_LuaModel_LIBRARY AND RBDL_LuaModel_INCLUDE_DIR)

IF (RBDL_URDFReader_LIBRARY AND RBDL_URDFReader_INCLUDE_DIR)
  SET (RBDL_URDFReader_FOUND TRUE)
ELSE(RBDL_URDFReader_LIBRARY AND RBDL_URDFReader_INCLUDE_DIR)
  MESSAGE(SEND_ERROR "Could not find RDBL ADDONS URDFReader")
ENDIF (RBDL_URDFReader_LIBRARY AND RBDL_URDFReader_INCLUDE_DIR)

IF (RBDL_CUSTOMFORCES_LIBRARY AND RBDL_CUSTOMFORCES_INCLUDE_DIR)
  SET (RBDL_GEOMETRY_FOUND TRUE)
ELSE (RBDL_CUSTOMFORCES_LIBRARY AND RBDL_CUSTOMFORCES_INCLUDE_DIR)
  MESSAGE(SEND_ERROR "Could not find RDBL ADDONS CUSTOM FORCES")
ENDIF (RBDL_CUSTOMFORCES_LIBRARY AND RBDL_CUSTOMFORCES_INCLUDE_DIR)

IF (RBDL_GEOMETRY_LIBRARY AND RBDL_GEOMETRY_INCLUDE_DIR)
  SET (RBDL_CUSTOMFORCES_FOUND TRUE)
ELSE (RBDL_GEOMETRY_LIBRARY AND RBDL_GEOMETRY_INCLUDE_DIR)
  MESSAGE(SEND_ERROR "Could not find RDBL ADDONS GEOMETRY")
ENDIF (RBDL_GEOMETRY_LIBRARY AND RBDL_GEOMETRY_INCLUDE_DIR)


IF (RBDL_FOUND)
   IF (NOT RBDL_FIND_QUIETLY)
      MESSAGE(STATUS "Found RBDL: ${RBDL_LIBRARY}")
   ENDIF (NOT RBDL_FIND_QUIETLY)

   foreach ( COMPONENT ${RBDL_FIND_COMPONENTS} )
     IF (RBDL_${COMPONENT}_FOUND)
       IF (NOT RBDL_FIND_QUIETLY)
         MESSAGE(STATUS "Found RBDL ${COMPONENT}: ${RBDL_${COMPONENT}_LIBRARY}")
       ENDIF (NOT RBDL_FIND_QUIETLY)
     ELSE (RBDL_${COMPONENT}_FOUND)
       MESSAGE(ERROR " Could not find RBDL ${COMPONENT}")
     ENDIF (RBDL_${COMPONENT}_FOUND)
   endforeach ( COMPONENT )
ENDIF (RBDL_FOUND)

MARK_AS_ADVANCED (
  RBDL_INCLUDE_DIR
  RBDL_LIBRARY
  RBDL_LuaModel_INCLUDE_DIR
  RBDL_LuaModel_LIBRARY
  RBDL_URDFReader_INCLUDE_DIR
  RBDL_URDFReader_LIBRARY
  RBDL_CUSTOMFORCES_INCLUDE_DIR
  RBDL_CUSTOMFORCES_LIBRARY
  RBDL_GEOMETRY_INCLUDE_DIR
  RBDL_GEOMETRY_LIBRARY
  )
