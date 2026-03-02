cmake_minimum_required(VERSION 3.16)
project(CATS_header_compile_test)

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

add_library(header_self_containment OBJECT
    ${HEADER_TEST_SOURCES}
)

target_include_directories(header_self_containment PRIVATE
    ${REPO_ROOT}/CATS
    ${REPO_ROOT}/CATS_Extentions
    ${REPO_ROOT}/DLM_CppTools
    ${REPO_ROOT}/DLM_MathTools
)

target_compile_features(header_self_containment PRIVATE cxx_std_11)
