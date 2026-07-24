cmake_minimum_required(VERSION 3.16)
project(UFBD CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)            # gnu++20

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(FBD_SRC ${CMAKE_SOURCE_DIR}/UFBD)

file(GLOB_RECURSE FBD_SOURCES CONFIGURE_DEPENDS ${FBD_SRC}/*.cpp)
list(FILTER FBD_SOURCES EXCLUDE REGEX "/Eigen/")

set(FBD_MAIN ${FBD_SRC}/main.cpp)
set(SBC_MAIN ${FBD_SRC}/Core/SBC/sbc_main.cpp)
set(FBD_LIB_SOURCES ${FBD_SOURCES})
list(REMOVE_ITEM FBD_LIB_SOURCES ${FBD_MAIN} ${SBC_MAIN})

add_library(fbdcore STATIC ${FBD_LIB_SOURCES})
target_include_directories(fbdcore PUBLIC
    ${FBD_SRC}
    ${FBD_SRC}/Models
    ${FBD_SRC}/Mcmc
    ${FBD_SRC}/Core/Parameters
    ${FBD_SRC}/Core/Tree
    ${FBD_SRC}/Core/IO
    ${FBD_SRC}/Core/Utils
    ${FBD_SRC}/Core/SBC
)

find_package(Threads REQUIRED)
target_link_libraries(fbdcore PUBLIC Threads::Threads)
target_compile_options(fbdcore PUBLIC $<$<NOT:$<CONFIG:Debug>>:-march=x86-64-v3>)

# -ffp-contract=off; ApproxLik. computeLnL compensated summation may break under contraction/fast-math
set_source_files_properties(${FBD_SRC}/Models/ApproxBranchLengthLikelihood.cpp
    PROPERTIES COMPILE_OPTIONS "-ffp-contract=off")

add_executable(ufbd ${FBD_MAIN})
target_link_libraries(ufbd PRIVATE fbdcore)
if(CMAKE_EXE_LINKER_FLAGS MATCHES "-static( |$)")
    set(_jem_suffixes "${CMAKE_FIND_LIBRARY_SUFFIXES}")
    set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
    find_library(JEMALLOC_STATIC NAMES jemalloc)
    set(CMAKE_FIND_LIBRARY_SUFFIXES "${_jem_suffixes}")
    set(JEMALLOC_LIB "${JEMALLOC_STATIC}")
else()
    find_library(JEMALLOC_SHARED NAMES jemalloc)
    set(JEMALLOC_LIB "${JEMALLOC_SHARED}")
endif()
if(JEMALLOC_LIB)
    target_link_libraries(ufbd PRIVATE ${JEMALLOC_LIB})
endif()

add_executable(sbc ${SBC_MAIN})
target_compile_definitions(sbc PRIVATE BUILD_SBC)
target_link_libraries(sbc PRIVATE fbdcore)
