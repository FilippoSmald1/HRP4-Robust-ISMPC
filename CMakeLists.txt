cmake_minimum_required(VERSION 2.8.12)
add_compile_options(-O3)
project(ISMCP)

find_package(DART 6.1.1 REQUIRED COMPONENTS utils-urdf gui-osg)

add_compile_options(-std=c++11)

include_directories(${DART_INCLUDE_DIRS})
include_directories(blasfeo/include)
include_directories(hpipm/include)

add_subdirectory(blasfeo)
add_subdirectory(hpipm)

file(GLOB srcs "*.cpp" "*.hpp" "*.h" "qpOASES/*.cpp" "*qpOASES/.hpp" "hpipm/include/*.h" "blasfeo/include/*.h")

add_executable(${PROJECT_NAME} ${srcs} ${APP_ICON_RESOURCE_WINDOWS})

target_link_libraries(${PROJECT_NAME} ${DART_LIBRARIES} blasfeo hpipm)
