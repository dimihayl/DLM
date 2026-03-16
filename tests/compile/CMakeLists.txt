cmake_minimum_required(VERSION 3.16)
project(CATS_header_compile_test)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
 set(CMAKE_CXX_FLAGS "-std=c++1y ${CMAKE_CXX_FLAGS}")
 set(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -undefined dynamic_lookup")
endif()

# Repo root
get_filename_component(REPO_ROOT
    ${CMAKE_CURRENT_SOURCE_DIR}/../..
    ABSOLUTE)

# Collect headers
file(GLOB_RECURSE PROJECT_HEADERS
    ${REPO_ROOT}/CATS/*.h
    ${REPO_ROOT}/CATS_Extentions/*.h
    ${REPO_ROOT}/DLM_CppTools/*.h
    ${REPO_ROOT}/DLM_MathTools/*.h
    ${REPO_ROOT}/CECA/*.h
    ${REPO_ROOT}/DLM_FemtoTools/*.h
    ${REPO_ROOT}/DLM_RootTools/*.h
    ${REPO_ROOT}/CATS_TUTORIAL_2019/*.h
)

set(HEADER_TEST_SOURCES)

foreach(header ${PROJECT_HEADERS})
    string(MD5 header_hash ${header})
    set(test_file ${CMAKE_BINARY_DIR}/header_${header_hash}.cpp)

    file(WRITE ${test_file}
        "#include \"${header}\"\n"
    )

    list(APPEND HEADER_TEST_SOURCES ${test_file})
endforeach()

SET(DLM_REPO "${CMAKE_CURRENT_LIST_DIR}/../..")
SET(CATS_INSTALL "${CMAKE_CURRENT_LIST_DIR}/../../install")

# Include user-provided paths
include("${CMAKE_SOURCE_DIR}/../../install/CMake/UserPaths.cmake")
include_directories(${GSL_INCLUDE})
link_directories(${GSL_LIB})
link_directories(${ROOT_PATH})

add_library(header_self_containment OBJECT
    ${HEADER_TEST_SOURCES}
)

target_include_directories(header_self_containment PRIVATE
    ${REPO_ROOT}/CATS
    ${REPO_ROOT}/CATS_Extentions
    ${REPO_ROOT}/DLM_CppTools
    ${REPO_ROOT}/DLM_MathTools
    ${REPO_ROOT}/CECA
    ${REPO_ROOT}/DLM_FemtoTools
    ${REPO_ROOT}/DLM_RootTools
    ${REPO_ROOT}/CATS_TUTORIAL_2019
)

execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --cflags" OUTPUT_VARIABLE CFLAGS)
execute_process(COMMAND bash -c "${ROOT_PATH}/bin/root-config --libs" OUTPUT_VARIABLE LIBS)
string(REGEX REPLACE "
$" "" CFLAGS "${CFLAGS}")
string(REGEX REPLACE "
$" "" LIBS "${LIBS}")
set(CFLAGS " -O2 -lgomp -pthread -fopenmp -std=c++11 -Wnonnull ${CFLAGS}")
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
set_target_properties( header_self_containment PROPERTIES COMPILE_FLAGS ${CFLAGS} ARCHIVE_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib" LIBRARY_OUTPUT_DIRECTORY "${CATS_INSTALL}/lib" RUNTIME_OUTPUT_DIRECTORY "${CATS_INSTALL}/bin")
target_link_libraries(  header_self_containment -L${GSL_LIB} -lgsl -lgslcblas ${LIBS} -lgomp)
target_compile_features(header_self_containment PRIVATE cxx_std_11)
