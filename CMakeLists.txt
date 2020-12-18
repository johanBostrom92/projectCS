cmake_minimum_required(VERSION 3.10.0)

# Specify project name
project(projectcs)

# Specify build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

# Make vim happy
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set extra compiler flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if(UNIX)
  set(CMAKE_CXX_FLAGS "-W -Wall")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif(UNIX)
if(WIN32)
  set(CMAKE_CXX_FLAGS "/EHsc")
endif(WIN32)

# Add source directories
aux_source_directory("${CMAKE_CURRENT_SOURCE_DIR}/src" PROJECT_SRCS)

# Add include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src")
file(GLOB_RECURSE PROJECT_HDRS src/*.hh)

# Specify library base directory
set(LIBRARY_DIR "${CMAKE_CURRENT_SOURCE_DIR}/lib")

# Define variable for linked libraries
set(PROJECT_LIBRARIES)

# Add executable for project
add_executable(${PROJECT_NAME} ${PROJECT_SRCS} ${PROJECT_HDRS})



#set(${PROJECT_LIBRARIES} ${PROJECT_LIBRARIES} ${PYTHON_LIBRARIES} )

# Link executable to libraries
target_link_libraries(${PROJECT_NAME} ${PROJECT_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()



# Build and use matplot++
add_subdirectory(lib/matplotplusplus-1.0.1)
target_link_libraries(${PROJECT_NAME} PUBLIC matplot)

# Link Python 
find_package(PythonLibs REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PUBLIC ${PYTHON_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dist")
