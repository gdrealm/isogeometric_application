set(CMAKE_INCLUDE_CURRENT_DIR ON)

set( GEOMOD_GEOMETRY_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/geomodcore.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/polygon.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/grid.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/tmesh/hilbertcurve.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/tmesh/simpletmesh.cpp
#    ${CMAKE_CURRENT_SOURCE_DIR}/geometry/tmesh.cpp
)

###############################################################
add_library(geomod STATIC ${GEOMOD_GEOMETRY_SRC} )
# target_link_libraries(superlu_mt ${BLAS_LIBRARIES} )
# install(TARGETS superlu DESTINATION libs )

