# Find Google Test
find_package(GTest REQUIRED)

# Add Homebrew lib and include directories
link_directories(/opt/homebrew/lib)
include_directories(/opt/homebrew/include)

# Enable RPATH on macOS
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_INSTALL_RPATH "/opt/homebrew/lib")
set(CMAKE_BUILD_WITH_INSTALL_RPATH ON)

# Define a function to add test targets with consistent settings
function(add_gtest_target TARGET_NAME SOURCE_FILE)
    add_executable(${TARGET_NAME} ${SOURCE_FILE})
    
    # Find the actual GTest libraries
    find_library(GTEST_LIBRARY_PATH gtest PATHS /opt/homebrew/lib NO_DEFAULT_PATH)
    find_library(GTEST_MAIN_LIBRARY_PATH gtest_main PATHS /opt/homebrew/lib NO_DEFAULT_PATH)
    
    if(GTEST_LIBRARY_PATH AND GTEST_MAIN_LIBRARY_PATH)
        message(STATUS "Found GTest libraries at ${GTEST_LIBRARY_PATH} and ${GTEST_MAIN_LIBRARY_PATH}")
    else()
        message(WARNING "GTest libraries not found in expected location")
    endif()
    
    # Link directly to the static libraries
    target_link_libraries(${TARGET_NAME}
        PRIVATE
        polygon_projection_core
        ${GTEST_LIBRARY_PATH}
        ${GTEST_MAIN_LIBRARY_PATH}
    )
    
    add_test(NAME ${TARGET_NAME} COMMAND ${TARGET_NAME})
endfunction()

# Add all test targets
add_gtest_target(polygon_tests polygon_tests.cpp)
add_gtest_target(projection_tests projection_tests.cpp)
add_gtest_target(polygon_comparison_tests polygon_comparison_tests.cpp)
add_gtest_target(analytical_comparison_tests analytical_comparison_tests.cpp) 
add_gtest_target(performance_benchmark_tests performance_benchmark_tests.cpp)
add_gtest_target(convexity_tests convexity_tests.cpp)