cmake_minimum_required(VERSION 3.12)

include(${CMAKE_CURRENT_LIST_DIR}/cmake/codipack-include.cmake)

project(
  ${CODIPACK_NAME}
  VERSION ${CODIPACK_VERSION}
  DESCRIPTION "CoDiPack - Code Differentiation Package"
  HOMEPAGE_URL "http://www.scicomp.uni-kl.de/software/codi")

add_library(${PROJECT_NAME} INTERFACE)

target_include_directories(
  ${PROJECT_NAME}
  INTERFACE $<BUILD_INTERFACE:${CODIPACK_INCLUDE_DIR}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

target_compile_features(${PROJECT_NAME} INTERFACE ${CODIPACK_CXX_VERSION})

include(GNUInstallDirs)
install(DIRECTORY ${CODIPACK_INCLUDE_DIR}/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY ${CODIPACK_CMAKE_DIR}/ DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${PROJECT_NAME}/cmake)
