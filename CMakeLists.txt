cmake_minimum_required(VERSION 3.0)

project(DEMO LANGUAGES CXX)

set (CMAKE_CXX_STANDARD 17)

include_directories(${PROJECT_SOURCE_DIR}/include)
link_directories(${PROJECT_SOURCE_DIR}/src)

add_library(${PROJECT_NAME} STATIC ${PROJECT_SOURCE_DIR}/src/problems.cpp)

target_include_directories(
    ${PROJECT_NAME}
    PUBLIC
        ${PROJECT_SOURCE_DIR}/src
)

foreach(LAB s2_lab2)
    set(LAB_ROOT_PATH src/)
    set(LAB_SOURCES
        ${LAB_ROOT_PATH}/problems.cpp
    )

    add_executable(main_dev_${LAB} ${LAB_ROOT_PATH}/main_dev.cpp ${LAB_SOURCES})
    
    target_include_directories(main_dev_${LAB} PUBLIC /opt/OpenBLAS/include)
    target_link_libraries(main_dev_${LAB} PRIVATE ${PROJECT_NAME} /opt/OpenBLAS/lib/libopenblas.a)
    target_compile_options(main_dev_${LAB} PUBLIC "-O3")

endforeach()
