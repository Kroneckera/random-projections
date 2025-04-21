cmake_minimum_required (VERSION 3.17)
project (random_projections VERSION 0.1.0 LANGUAGES CXX)

# ──────────────────────────  Global compiler settings  ──────────────────────────
set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_POSITION_INDEPENDENT_CODE ON)

# ───────────────────────────────  Dependencies  ────────────────────────────────
find_package (pybind11 CONFIG REQUIRED)          # provides pybind11::headers, …
# pybind11::headers is header‑only → no link libs pulled into the core archive

# ───────────────────────────────  Paths / names  ───────────────────────────────
set (PY_PACKAGE      projection)
set (PY_PACKAGE_DIR  ${CMAKE_CURRENT_SOURCE_DIR}/python/${PY_PACKAGE})

# ─────────────────────────  Core static C++ library  ───────────────────────────
file (GLOB_RECURSE CPP_SOURCES CONFIGURE_DEPENDS  cpp/src/*.cpp)

# Drop every file that contains bindings or is flagged “minimal”
list (FILTER CPP_SOURCES EXCLUDE REGEX "python_bindings\\.cpp$")
list (FILTER CPP_SOURCES EXCLUDE REGEX "minimal_.*\\.cpp$")

add_library (random_projections_core STATIC ${CPP_SOURCES})
target_include_directories (random_projections_core PUBLIC cpp/include)
target_link_libraries      (random_projections_core PUBLIC pybind11::headers)

# ────────────────────────  Python extension module  ────────────────────────────
# Builds _core.(so|pyd|dll) from the single bindings file
pybind11_add_module (_core cpp/src/python_bindings.cpp)
target_link_libraries (_core PRIVATE random_projections_core)

set_target_properties (_core PROPERTIES
    PREFIX ""                                  # remove "lib" → _core.so / .pyd
    LIBRARY_OUTPUT_DIRECTORY  ${PY_PACKAGE_DIR}
    RUNTIME_OUTPUT_DIRECTORY  ${PY_PACKAGE_DIR}   # Windows
    ARCHIVE_OUTPUT_DIRECTORY  ${PY_PACKAGE_DIR})  # MSVC import libs

# On MSVC the correct extension for Python modules is .pyd
if (MSVC)
    set_target_properties (_core PROPERTIES SUFFIX ".pyd")
endif()

# ───────────────────────────────  Installation  ────────────────────────────────
install (TARGETS _core DESTINATION ${PY_PACKAGE_DIR})
