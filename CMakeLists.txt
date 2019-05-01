cmake_minimum_required(VERSION 2.8)
project(SimpleDenoise LANGUAGES C)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
SET(CMAKE_BUILD_TYPE "Release")

include_directories(include)
add_executable(${PROJECT_NAME} main.c)

target_link_libraries(${PROJECT_NAME} -lm)