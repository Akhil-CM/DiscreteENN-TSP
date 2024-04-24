cmake_minimum_required(VERSION 3.12)
# set the project name
project(DiscreteENN_TSP VERSION 0.1)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_C_COMPILER "clang")
set(CMAKE_CXX_COMPILER "clang++")

# # we default to Release build type
# if(NOT CMAKE_BUILD_TYPE)
#     set(CMAKE_BUILD_TYPE "Debug")
# endif()

if (MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(-Wall -Wextra -pedantic)
endif()

include_directories( "${PROJECT_SOURCE_DIR}/include" )

# add_library(utilities OBJECT "${PROJECT_SOURCE_DIR}/src/utilities.cxx")

# add the executable
add_executable(enn_tsp "${PROJECT_SOURCE_DIR}/main.cpp" "${PROJECT_SOURCE_DIR}/tsp_discrete_enn.hpp" "${PROJECT_SOURCE_DIR}/utils.hpp")

# Set compiler options for different build types
target_compile_options(enn_tsp PRIVATE
    $<$<CONFIG:Release>:-O3;-march=native;-mtune=native>
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,float-divide-by-zero,nullability;-fno-omit-frame-pointer;-Og;-g>
)

# Link options, needed specifically for the sanitizer in Debug mode
target_link_options(enn_tsp PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,float-divide-by-zero,nullability;-fno-omit-frame-pointer>
)

target_link_libraries(enn_tsp PRIVATE
		      sfml-graphics
                      sfml-window
                      sfml-system
		      )

install(TARGETS enn_tsp DESTINATION "${PROJECT_SOURCE_DIR}/bin")

# Copy to source directory
add_custom_target(
    copy-compile-commands ALL
    DEPENDS
    ${CMAKE_SOURCE_DIR}/compile_commands.json
    )
add_custom_command(
    OUTPUT ${CMAKE_SOURCE_DIR}/compile_commands.json
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${CMAKE_BINARY_DIR}/compile_commands.json
        ${CMAKE_SOURCE_DIR}/compile_commands.json
    DEPENDS
        # Unlike "proper" targets like executables and libraries,
        # custom command / target pairs will not set up source
        # file dependencies, so we need to list file explicitly here
        generate-compile-commands
        ${CMAKE_BINARY_DIR}/compile_commands.json
    )

# Generate the compilation commands. Necessary so cmake knows where it came
# from and if for some reason you delete it.
add_custom_target(generate-compile-commands
    DEPENDS
        ${CMAKE_BINARY_DIR}/compile_commands.json
    )
add_custom_command(
    OUTPUT ${CMAKE_BINARY_DIR}/compile_commands.json
    COMMAND ${CMAKE_COMMAND} -B${CMAKE_BINARY_DIR} -S${CMAKE_SOURCE_DIR}
    )
