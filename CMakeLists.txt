cmake_minimum_required(VERSION 3.1.0)
project(complex_simplification)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

if(MSVC)
  # Enable parallel compilation
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP /bigobj")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR} )
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR} )
else()
  # Libigl requires a modern C++ compiler that supports c++11
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "." )
endif()

option(LIBIGL_USE_STATIC_LIBRARY "Use LibIGL as static library" ON)

# add a customizable menu bar
option(LIBIGL_WITH_NANOGUI     "Use Nanogui menu"   ON)

# libigl options: choose your dependencies (by default everything is OFF except opengl)
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_GLFW        "Use GLFW"           ON)
option(LIBIGL_WITH_BBW         "Use BBW"            OFF)
option(LIBIGL_WITH_EMBREE      "Use Embree"         ON)
option(LIBIGL_WITH_PNG         "Use PNG"            ON)
option(LIBIGL_WITH_TETGEN      "Use Tetgen"         OFF)
option(LIBIGL_WITH_TRIANGLE    "Use Triangle"       OFF)
option(LIBIGL_WITH_XML         "Use XML"            OFF)
option(LIBIGL_WITH_LIM         "Use LIM"            OFF)
option(LIBIGL_WITH_COMISO      "Use CoMiso"         OFF)
option(LIBIGL_WITH_MATLAB      "Use Matlab"         OFF) # This option is not supported yet
option(LIBIGL_WITH_MOSEK       "Use MOSEK"          OFF) # This option is not supported yet
option(LIBIGL_WITH_CGAL        "Use CGAL"           OFF)
if(LIBIGL_WITH_CGAL) # Do not remove or move this block, the cgal build system fails without it
  find_package(CGAL REQUIRED)
  set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE CACHE BOOL "CGAL's CMAKE Setup is super annoying ")
  include(${CGAL_USE_FILE})
endif()

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl/shared/cmake" "libigl")

# libigl information
message("libigl includes: ${LIBIGL_INCLUDE_DIRS}")
message("libigl libraries: ${LIBIGL_LIBRARIES}")
message("libigl extra sources: ${LIBIGL_EXTRA_SOURCES}")
message("libigl extra libraries: ${LIBIGL_EXTRA_LIBRARIES}")
message("libigl definitions: ${LIBIGL_DEFINITIONS}")

# Build TBB
set(TBB_BUILD_STATIC ON CACHE BOOL " " FORCE)
set(TBB_BUILD_SHARED OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TBBMALLOC_PROXY OFF CACHE BOOL " " FORCE)
set(TBB_BUILD_TESTS OFF CACHE BOOL " " FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/tbb extern_build/tbb)
set_property(TARGET tbb_static tbb_def_files PROPERTY FOLDER "dependencies")

# Build VCG
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/vcg extern_build/vcg)
set_property(TARGET vcg PROPERTY FOLDER "dependencies")
# Prepare the build environment
include_directories(
${LIBIGL_INCLUDE_DIRS}
  # Intel Thread Building Blocks
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/tbb/include
  # VCG
  ${CMAKE_CURRENT_SOURCE_DIR}/extern/vcg
  )
add_definitions(${LIBIGL_DEFINITIONS})
add_definitions(-DENABLE_SERIALIZATION)


file(GLOB header *.h)
file(GLOB source *.cpp)

add_executable(complex_simplification ${source} ${header})
target_link_libraries(${PROJECT_NAME}  tbb_static vcg ${LIBIGL_LIBRARIES} ${LIBIGL_EXTRA_LIBRARIES})