cmake_minimum_required(VERSION 3.0)

project(csv_analyzer LANGUAGES C)

add_executable(csv_analyzer main.c)

install(TARGETS csv_analyzer RUNTIME DESTINATION bin)
