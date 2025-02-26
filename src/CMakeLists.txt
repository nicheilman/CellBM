cmake_minimum_required(VERSION 3.0)

#Add the created libraries in the project

if(ENABLE_PYTHON_BINDINGS)
    add_subdirectory(python_bindings)
endif()


FILE(GLOB SRC "*.cpp")

#Create the source library
add_library(src ${SRC})

#Add all the header files required to run the root file of the project
target_include_directories(src PUBLIC 
    "${PROJECT_SOURCE_DIR}/include"
)

#Locate the openMP library to parallelize the code
find_package(OpenMP REQUIRED)

#Link root library of the src directory with all the other libraries created
target_link_libraries(src PUBLIC    
    OpenMP::OpenMP_CXX          
)


