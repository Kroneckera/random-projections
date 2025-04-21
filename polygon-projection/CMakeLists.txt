cmake_minimum_required(VERSION 3.14)
project(polygon_projection VERSION 0.1.0 LANGUAGES CXX)

# C++ standard settings - MUST use C++17 for pybind11
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Add compiler warnings based on compiler
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-Wall -Wextra -Wpedantic)
elseif(MSVC)
    add_compile_options(/W4)
endif()

# Platform-specific settings
if(APPLE)
    add_compile_options(-std=c++17)
elseif(UNIX AND NOT APPLE)
    add_compile_options(-std=c++17)
elseif(WIN32)
    # Disable warnings about standard C functions on Windows
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
    # Make sure static libraries are generated with correct flags
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
    # Don't warn about standard-conforming behavior
    if(MSVC)
        add_compile_options(/permissive-)
    endif()
endif()

# Include directories for headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add the core library
add_library(polygon_projection_core
    src/polygon.cpp
    src/projection.cpp
    src/interior_projection.cpp
    src/boundary_projection.cpp
    src/sampling.cpp
    src/integration.cpp
    src/utils.cpp
)

# Set properties on the core library
target_include_directories(polygon_projection_core PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Make sure the core library is built with position-independent code
set_target_properties(polygon_projection_core PROPERTIES
    POSITION_INDEPENDENT_CODE ON
)

# Python bindings (optional)
option(BUILD_PYTHON_BINDINGS "Build Python bindings" ON)
if(BUILD_PYTHON_BINDINGS)
    message(STATUS "Python bindings enabled")
    
    # Try to find Python
    if(DEFINED PYTHON_EXECUTABLE)
        message(STATUS "Using specified Python executable: ${PYTHON_EXECUTABLE}")
    else()
        # Find Python interpreter
        find_program(PYTHON_EXECUTABLE 
            NAMES python3 python
            DOC "Python interpreter"
        )
        
        if(NOT PYTHON_EXECUTABLE)
            message(FATAL_ERROR "Python executable not found. Specify with -DPYTHON_EXECUTABLE=/path/to/python")
        endif()
        
        message(STATUS "Using detected Python executable: ${PYTHON_EXECUTABLE}")
    endif()
    
    # Get Python version
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
        OUTPUT_VARIABLE PYTHON_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    message(STATUS "Building for Python version: ${PYTHON_VERSION}")
    
    # Find Python development libraries
    find_package(Python ${PYTHON_VERSION} COMPONENTS Interpreter Development REQUIRED)
    
    message(STATUS "Python includes: ${Python_INCLUDE_DIRS}")
    message(STATUS "Python libraries: ${Python_LIBRARIES}")
    
    # Check for pybind11
    execute_process(
        COMMAND ${PYTHON_EXECUTABLE} -c "import pybind11; print(pybind11.get_include())"
        RESULT_VARIABLE PYBIND11_IMPORT_RESULT
        OUTPUT_VARIABLE PYBIND11_INCLUDE_DIR
        ERROR_QUIET
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    
    if(NOT PYBIND11_IMPORT_RESULT EQUAL 0)
        message(STATUS "pybind11 not found, attempting to install it...")
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -m pip install pybind11
            RESULT_VARIABLE PIP_RESULT
        )
        
        if(NOT PIP_RESULT EQUAL 0)
            message(FATAL_ERROR "Failed to install pybind11. Please install it manually with: pip install pybind11")
        endif()
        
        # Try again to get the include path
        execute_process(
            COMMAND ${PYTHON_EXECUTABLE} -c "import pybind11; print(pybind11.get_include())"
            RESULT_VARIABLE PYBIND11_IMPORT_RESULT
            OUTPUT_VARIABLE PYBIND11_INCLUDE_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()
    
    # Check if we found pybind11
    if(PYBIND11_INCLUDE_DIR)
        message(STATUS "Found pybind11: ${PYBIND11_INCLUDE_DIR}")
        
        # Set the Python module output directory
        set(PYTHON_MODULE_OUTPUT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/python/polygon_projection")
        
        # Make sure the output directory exists
        file(MAKE_DIRECTORY ${PYTHON_MODULE_OUTPUT_DIR})
        
        # Create the Python extension module
        add_library(_core MODULE src/python_bindings.cpp)
        target_link_libraries(_core PRIVATE polygon_projection_core)
        target_include_directories(_core PRIVATE ${PYBIND11_INCLUDE_DIR} ${Python_INCLUDE_DIRS})
        
        # Platform-specific settings for the extension module
        if(WIN32)
            set_target_properties(_core PROPERTIES SUFFIX ".pyd")
        else()
            set_target_properties(_core PROPERTIES SUFFIX ".so")
        endif()
        
        # Always link against Python libraries
        target_link_libraries(_core PRIVATE ${Python_LIBRARIES})
        
        # Set the output name and location
        set_target_properties(_core PROPERTIES
            PREFIX ""
            OUTPUT_NAME "_core"
            LIBRARY_OUTPUT_DIRECTORY "${PYTHON_MODULE_OUTPUT_DIR}"
        )
        
        message(STATUS "Python module will be output to: ${PYTHON_MODULE_OUTPUT_DIR}")
    else()
        message(FATAL_ERROR "pybind11 not found. Python bindings cannot be built")
    endif()
endif()

# Testing (optional)
option(BUILD_TESTING "Build tests" ON)
if(BUILD_TESTING)
    enable_testing()
    
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/CMakeLists.txt")
        add_subdirectory(tests/cpp)
        message(STATUS "C++ tests will be built")
    else()
        message(STATUS "C++ tests directory not found, skipping tests")
    endif()
endif()

# Examples (optional)
option(BUILD_EXAMPLES "Build example applications" ON)
if(BUILD_EXAMPLES)
    # Check for each example file and build it if it exists
    set(EXAMPLE_FILES
        "basic_example.cpp"
        "simple_full_example.cpp"
    )
    
    foreach(EXAMPLE_FILE ${EXAMPLE_FILES})
        if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples/${EXAMPLE_FILE}")
            get_filename_component(EXAMPLE_NAME ${EXAMPLE_FILE} NAME_WE)
            add_executable(${EXAMPLE_NAME} examples/${EXAMPLE_FILE})
            target_link_libraries(${EXAMPLE_NAME} PRIVATE polygon_projection_core)
            message(STATUS "Example ${EXAMPLE_NAME} will be built")
        endif()
    endforeach()
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "=== ${PROJECT_NAME} Configuration Summary ===")
message(STATUS "Version:           ${PROJECT_VERSION}")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Python bindings:   ${BUILD_PYTHON_BINDINGS}")
message(STATUS "Tests:             ${BUILD_TESTING}")
message(STATUS "Examples:          ${BUILD_EXAMPLES}")
message(STATUS "=====================================")
message(STATUS "")