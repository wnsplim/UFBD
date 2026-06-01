cmake_minimum_required(VERSION 3.16)
project(FossilizedBirthDeathProcesses CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS ON)            # gnu++20

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

set(FBD_SRC ${CMAKE_SOURCE_DIR}/FossilizedBirthDeathProcesses)

file(GLOB_RECURSE FBD_SOURCES CONFIGURE_DEPENDS ${FBD_SRC}/*.cpp)
list(FILTER FBD_SOURCES EXCLUDE REGEX "/Eigen/")

add_executable(fbd ${FBD_SOURCES})

target_include_directories(fbd PRIVATE
    ${FBD_SRC}                          # bundled Eigen ("Eigen/Dense") + root-level headers
    ${FBD_SRC}/Models
    ${FBD_SRC}/Mcmc
    ${FBD_SRC}/Core/Parameters
    ${FBD_SRC}/Core/Tree
    ${FBD_SRC}/Core/IO
    ${FBD_SRC}/Core/Utils
)

find_package(Threads REQUIRED)
target_link_libraries(fbd PRIVATE Threads::Threads)
