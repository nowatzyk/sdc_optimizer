cmake_minimum_required(VERSION 3.10)
project(csv_analyzer LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(FLEX REQUIRED)
find_package(BISON REQUIRED)

BISON_TARGET(
    pragma_parser
    parser.y
    ${CMAKE_CURRENT_BINARY_DIR}/parser.c
    DEFINES_FILE ${CMAKE_CURRENT_BINARY_DIR}/parser.h
)

FLEX_TARGET(
    pragma_lexer
    lexer.l
    ${CMAKE_CURRENT_BINARY_DIR}/lexer.c
    COMPILE_FLAGS "--header-file=${CMAKE_CURRENT_BINARY_DIR}/lex.yy.h"
)

ADD_FLEX_BISON_DEPENDENCY(pragma_lexer pragma_parser)

add_executable(csv_analyzer
    main.cpp
    lsq_fit.c
    fit_functions.c
    xrand.cpp
    parameter.cpp
    eval_cache.cpp
    expression.cpp sum_expression.cpp geomean_expression.cpp min_expression.cpp max_expression.cpp
    loop_complex.cpp lsq_fit_function.cpp anneal.cpp
    nodes_of_interest.cpp spice_deck.cpp parser_interface.cpp
    ${BISON_pragma_parser_OUTPUTS}
    ${FLEX_pragma_lexer_OUTPUTS}
)

target_include_directories(csv_analyzer PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS csv_analyzer RUNTIME DESTINATION bin)
