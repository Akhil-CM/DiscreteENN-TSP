cmake_minimum_required(VERSION 3.15)
# set the project name
project(DiscreteENN_TSP VERSION 1.0
    DESCRIPTION
    "A fast algorithm to find reasonable solutions to TSP"
    LANGUAGES CXX)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with "
                        "CMakeLists.txt file). Please make a build subdirectory. Feel free to "
                        "remove CMakeCache.txt and CMakeFiles.")
endif()

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

file(GLOB HEADER_LIST CONFIGURE_DEPENDS
    "${PROJECT_SOURCE_DIR}/include/*.hpp")

add_library(tsp_enn_lib)
target_sources(tsp_enn_lib
  PRIVATE
    "${PROJECT_SOURCE_DIR}/src/tsp_discrete_enn.cxx"
  PUBLIC
    "${HEADER_LIST}"
)

find_program(ROOT_CONFIG_EXEC root-config)

if(NOT ROOT_CONFIG_EXEC)
    message(FATAL_ERROR "root-config not found")
endif()

execute_process(
    COMMAND ${ROOT_CONFIG_EXEC} --cflags
    OUTPUT_VARIABLE ROOT_CXX_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND ${ROOT_CONFIG_EXEC} --libs
    OUTPUT_VARIABLE ROOT_LD_FLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Separate the flags into a CMake list
separate_arguments(ROOT_CXX_FLAGS_LIST UNIX_COMMAND "${ROOT_CXX_FLAGS}")
separate_arguments(ROOT_LD_FLAGS_LIST UNIX_COMMAND "${ROOT_LD_FLAGS}")

# Set the target properties
target_compile_options(tsp_enn_lib PRIVATE ${ROOT_CXX_FLAGS_LIST})
target_link_libraries(tsp_enn_lib PRIVATE ${ROOT_LD_FLAGS_LIST})

target_include_directories(tsp_enn_lib PUBLIC "${PROJECT_SOURCE_DIR}/include")
target_link_libraries(tsp_enn_lib PRIVATE
		      sfml-graphics
                      sfml-window
                      sfml-system
		      )

add_executable(tsp_enn "${PROJECT_SOURCE_DIR}/src/main.cpp" "${HEADER_LIST}")
target_include_directories(tsp_enn PUBLIC "${PROJECT_SOURCE_DIR}/include")

# target_compile_features(enn_tsp PRIVATE cxx_std_17)
target_compile_options(tsp_enn PRIVATE
    $<$<CONFIG:Release>:-O3;-march=native;-mtune=native>
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,float-divide-by-zero,nullability;-fno-omit-frame-pointer;-Og;-g>
)
target_link_options(tsp_enn PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address,undefined,float-divide-by-zero,nullability;-fno-omit-frame-pointer>
)

target_link_libraries(tsp_enn PRIVATE tsp_enn_lib)

install(TARGETS tsp_enn DESTINATION "${PROJECT_SOURCE_DIR}/bin")

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

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    add_custom_target(generate_git_commit ALL
        DEPENDS "${CMAKE_BINARY_DIR}/last_git_commit.txt")

    add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/last_git_commit.txt"
        COMMAND ${GIT_EXECUTABLE} log -n 2 > "${CMAKE_BINARY_DIR}/last_git_commit.txt" 2>&1)
endif()

# find_program(BASH bash HINTS /bin)
# find_package(Git QUIET)
# if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
#     add_custom_target(generate_git_commit ALL
#         DEPENDS "${CMAKE_BINARY_DIR}/last_git_commit.txt")

#     add_custom_command(OUTPUT "${CMAKE_BINARY_DIR}/last_git_commit.txt"
#         COMMAND ${CMAKE_COMMAND} -E env ${BASH} -c "${GIT_EXECUTABLE} log -n 2 > ${CMAKE_BINARY_DIR}/last_git_commit.txt 2>&1")
#         VERBATIM)
# endif()
