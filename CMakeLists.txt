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

add_executable(ufbd ${FBD_MAIN})
target_link_libraries(ufbd PRIVATE fbdcore)
find_library(JEMALLOC_LIB NAMES jemalloc)
if(JEMALLOC_LIB)
    target_link_libraries(ufbd PRIVATE ${JEMALLOC_LIB})
endif()

add_executable(sbc ${SBC_MAIN})
target_compile_definitions(sbc PRIVATE BUILD_SBC)
target_link_libraries(sbc PRIVATE fbdcore)
