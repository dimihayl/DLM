CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#It is recommended that you copy all CMake related files in a local folder of yours,
#else whenever you pull/push there will be conflicts to resolve

project(PION_FITTER)
# SET PATHS #
SET(PROJECT_DESTINATION "/home/dmihaylov/Temp/CATS_FITTER1")#the destination of your project
SET(ROOT_PATH "/home/dmihaylov/Apps/root-6.14.00/obj")#path to ROOT
#SET(ROOT_PATH "/home/dmihaylov/root")
SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
SET(GSL_LIB "/usr/lib")#where are the GSL .a and .so files
SET(CATS_PATH "/home/dmihaylov/Apps/CATS/bin")#the CATS main folder (containing cats-config)
SET(CATS_TYPE "dev")#basic, extended or dev

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DESTINATION})
set_target_properties(PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_DESTINATION}/bin")

add_executable(PION_FITTER main.cpp
CATS_ExampleFitter1.cpp
#ADD BELOW ANY OTHER .CPP FILES THAT YOU WOULD LIKE TO COMPILE (FULL PATH)
)

# DO NOT CHANGE THE REST #

execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --cflags" OUTPUT_VARIABLE CFLAGS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --libs" OUTPUT_VARIABLE LIBS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --glibs" OUTPUT_VARIABLE GLIBS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --incdir" OUTPUT_VARIABLE ROOT_INCLUDE)
execute_process(COMMAND bash -c "${CATS_PATH}/cats-config --incdir" OUTPUT_VARIABLE CATS_INCLUDE)
execute_process(COMMAND bash -c "${CATS_PATH}/cats-config --libs-${CATS_TYPE}" OUTPUT_VARIABLE CATS_LIBS)

string(REGEX REPLACE "\n$" "" CFLAGS "${CFLAGS}")
string(REGEX REPLACE "\n$" "" LIBS "${LIBS}")
string(REGEX REPLACE "\n$" "" GLIBS "${GLIBS}")
string(REGEX REPLACE "\n$" "" ROOT_INCLUDE "${ROOT_INCLUDE}")
string(REGEX REPLACE "\n$" "" CATS_INCLUDE "${CATS_INCLUDE}")
string(REGEX REPLACE "\n$" "" CATS_LIBS "${CATS_LIBS}")

string(APPEND CFLAGS " -O3")

set_target_properties(PION_FITTER PROPERTIES COMPILE_FLAGS ${CFLAGS})

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

message(STATUS ${PROJECT_SOURCE_DIR})
message(STATUS ${CATS_INCLUDE})
message(STATUS ${GSL_INCLUDE})
message(STATUS ${ROOT_INCLUDE})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${ROOT_INCLUDE})

target_link_libraries(PION_FITTER -L${CATS_LIBS} ${LIBS} -L${GSL_LIB} -lgsl -lgslcblas)

