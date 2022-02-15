cmake_minimum_required(VERSION 3.10)

project(rtg_library) #Robot Trajectory Generation (RTG) Library 
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -pedantic")

find_package(Eigen3 REQUIRED)
   
include_directories(
	${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${EIGEN3_INCLUDE_DIR}
)

file(GLOB SRC_LIST
        "${PROJECT_SOURCE_DIR}/include/*.h"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        )
include(${CMAKE_CURRENT_SOURCE_DIR}/visulization/CMakeLists.txt)

#add_executable(polynomial ${SRC_LIST} main.cpp)

# Executable
add_executable(polynomial ${SRC_LIST} )
target_link_libraries(polynomial PRIVATE matplotlib_cpp)
set_target_properties(polynomial PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")




