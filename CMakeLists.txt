cmake_minimum_required(VERSION 3.0)

project(csv_analyzer)

add_executable(csv_analyzer main.cpp)

install(TARGETS csv_analyzer RUNTIME DESTINATION bin)
