CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

project(PION_FITTER)
# SET PATHS #
SET(PROJECT_DESTINATION "/home/dmihaylov/Dudek_Ubuntu/Temp/cat/ProjectSource")#the destination of your project
SET(ROOT_INCLUDE "/home/dmihaylov/Apps/root-6.14.00/obj/include")#where are all ROOT related .h files
SET(ROOT_LIB "/home/dmihaylov/Apps/root-6.14.00/obj/lib")#where are the ROOT related .a and .so files
SET(GSL_INCLUDE "/usr/include/gsl")#where are all GSL related .h files
SET(GSL_LIB "/usr/lib")#where are the GSL .a and .so files
SET(CATS_INCLUDE "/home/dmihaylov/Dudek_Ubuntu/Temp/cat/CatsSource/include")#where are all CATS related .h files
SET(CATS_LIB "/home/dmihaylov/Dudek_Ubuntu/Temp/cat/CatsSource/lib")#where are the CATS related .a files

SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_DESTINATION})
set_target_properties(PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_DESTINATION}/bin")

add_executable(PION_FITTER main.cpp
CATS_ExampleFitter1.cpp
#ADD BELOW ANY OTHER .CPP FILES THAT YOU WOULD LIKE TO COMPILE (FULL PATH)
)
#N.B. Depending on your ROOT version, you might sill want to change the target_link_libraries at the end of the file
# DO NOT CHANGE THE REST #

set_target_properties(PION_FITTER PROPERTIES COMPILE_FLAGS "-pthread -std=c++11 -m64")

SET(VERSION_MAJOR 1)
SET(VERSION_MINOR 0)
SET(VERSION_PATCH 0)
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
include_directories(${CATS_INCLUDE})
include_directories(${GSL_INCLUDE})
include_directories(${ROOT_INCLUDE})

#use this for ROOT5
#target_link_libraries(PION_FITTER -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${ROOT_LIB} -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -L${GSL_LIB} -lgsl -lgslcblas)

#use this for ROOT6
target_link_libraries(PION_FITTER -L${CATS_LIB} -lCATSbasic -lCATSextended -lCATSdev -L${ROOT_LIB} -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -L${GSL_LIB} -lgsl -lgslcblas)


