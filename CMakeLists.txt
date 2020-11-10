cmake_minimum_required(VERSION 3.10.0)

# Specify project name
project(projectcs)

# Specify build type
set(CMAKE_BUILD_TYPE Debug)

# Make vim happy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set extra compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(UNIX)
  set(CMAKE_CXX_FLAGS "-W -Wall")
endif(UNIX)
if(WIN32)
  set(CMAKE_CXX_FLAGS "/EHsc")
endif(WIN32)

# Add source directories
aux_source_directory("${CMAKE_CURRENT_SOURCE_DIR}/src" PROJECT_SRCS)

# Add include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")

# Specify library base directory
set(LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Define variable for linked libraries
set(PROJECT_LIBRARIES)

# Add executable for project
add_executable(${PROJECT_NAME} ${PROJECT_SRCS})

# Link executable to libraries
target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()

set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dist")