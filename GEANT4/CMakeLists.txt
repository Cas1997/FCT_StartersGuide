cmake_minimum_required(VERSION 3.16...3.27)

project(FCTV2)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

file(GLOB MACRO_FILES "*.mac")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

file(GLOB MACRO_FILES "*.cfg")
file(COPY ${MACRO_FILES} DESTINATION ${PROJECT_BINARY_DIR})

add_executable(FCT run.cc ${sources} ${headers})
target_link_libraries(FCT ${Geant4_LIBRARIES})

add_custom_target(FCTV2 DEPENDS FCT)