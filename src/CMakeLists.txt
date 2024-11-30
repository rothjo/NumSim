cmake_minimum_required(VERSION 3.8)

# Define the project name.
project(numsim)

# Specify the name of the executable (${PROJECT_NAME} which is equal to what was set in the project() command).
# Also specify the source files.
add_executable(${PROJECT_NAME}
   main.cpp
   storage/array2d.cpp
   storage/FieldVariable.cpp
   discretization/centralDifferences.cpp
   discretization/discretization.cpp
   discretization/staggeredGrid.cpp
   discretization/donorCell.cpp
   pressureSolver/pressureSolver.cpp
   pressureSolver/gaussSeidel.cpp
   pressureSolver/sor.cpp
   pressureSolver/parallelGaussSeidel.cpp
   pressureSolver/parallel_sor.cpp
   pressureSolver/parallelPressureSolver.cpp
   pressureSolver/cg.cpp
   computation/computation.cpp
   computation/parallelComputation.cpp
   settings.cpp
   output_writer/output_writer.cpp
   output_writer/output_writer_text.cpp
   output_writer/output_writer_paraview.cpp
   output_writer/output_writer_paraview_parallel.cpp
   output_writer/output_writer_text_parallel.cpp
   partitioning/partitioning.cpp

)

# Add the project directory to include directories, to be able to include all project header files from anywhere
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR})

# Search for the external package "VTK"
# where is vtk? ask chatgpt
# set(VTK_DIR /home/jonas/software/VTK-7.1.1/build)
find_package(VTK)

# Output various cmake variables for demonstration purpose
message("If VTK was found on the system: VTK_FOUND: ${VTK_FOUND}")
message("The directory of VTK: VTK_DIR: ${VTK_DIR}")
message("The include directory of VTK: VTK_INCLUDE_DIRS: ${VTK_INCLUDE_DIRS}")

# If an installation of vtk was found on the system
if (VTK_FOUND)
  include_directories(${VTK_INCLUDE_DIRS})               # add the include directory where the header files are for the compiler
  target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES}) # add the libraries for the linker
endif(VTK_FOUND)

find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

find_package(Java REQUIRED)
find_package(JNI REQUIRED)

include_directories(${JNI_INCLUDE_DIRS})
link_directories(/usr/lib/jvm/java-11-openjdk-amd64/lib/)

target_link_libraries(${PROJECT_NAME} ${JNI_LIBRARIES})


#export JAVA_HOME=/usr/lib/jvm/java-11-openjdk-amd64
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-11-openjdk-amd64/lib/server
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/jvm/java-11-openjdk-amd64/lib/

# install numsim executable in build directory
install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/../build)

# Add additonial compile options to enable more warnings
add_compile_options(-Wall -Wextra)

# Set the version of the C++ standard to use, we use C++14, published in 2014
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message("CMAKE_BUILD_TYPE: ${CMAKE_BUILD_TYPE}")