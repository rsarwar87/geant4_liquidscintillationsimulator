#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(LiquidScintillatorSimulator)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(LiquidScintillatorSimulator LiquidScintillatorSimulator.cc ${sources} ${headers})
target_link_libraries(LiquidScintillatorSimulator ${Geant4_LIBRARIES} ${ARMADILLO_LIBRARIES})
target_link_libraries(LiquidScintillatorSimulator /media/2TB/workspace/rsarwar/server/fission_v2.0.5/lib/libFission.so)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build LiquidScintillatorSimulator. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(LiquidScintillatorSimulator_SCRIPTS
    LiquidScintillatorSimulator.out
    LiquidScintillatorSimulator.in
    optPhoton.mac
    gui.mac
    icons.mac
    run.png
    vis.mac
  )

foreach(_script ${LiquidScintillatorSimulator_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS LiquidScintillatorSimulator DESTINATION bin)

