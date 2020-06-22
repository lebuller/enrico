cmake_minimum_required (VERSION 3.8)

# Set cxx standard and compiler flags
set( CMAKE_CXX_STANDARD 11)
set( CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Xlinker --no-as-needed -Xlinker --add-needed")

# Sets link to the directory where the foam lib*.so files are stored
link_directories($ENV{FOAM_LIBBIN})

# Adds compiler flags needed in OpenFOAM
add_definitions( -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository)


# Add source to this project's executable
# Mimics the options in Make/files file
add_library (openFOAM-imported enricoFoamLibrary.C solid/solidRegionDiffNo.C fluid/compressibleCourantNo.C)

#Adds the include directories
# Mimics the options in the Make/options file
include_directories(
    ./
    ./fluid/
    ./solid/
    ./porousFluid/
    ./porousSolid/
    ./include/
    $ENV{WM_PROJECT_DIR}/src/meshTools/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/sampling/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/finiteVolume/cfdTools/
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/basic/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/specie/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/reactionThermo/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/solidThermo/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/ODE/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/combustionModels/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/momentumTransportModels/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/MomentumTransportModels/compressible/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/ThermophysicalTransportModels/rhoReactionThermo/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/radiationModels/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/regionModels/regionModel/lnInclude/

    $ENV{WM_PROJECT_DIR}/src/OpenFOAM/lnInclude/
    $ENV{WM_PROJECT_DIR}/src/OSspecific/POSIX/lnInclude/
)

set( EXTRA_LIBS dl m )
set( LIBS
          fluidThermophysicalModels
          specie
          reactionThermophysicalModels
          solidThermo
          chemistryModel
          ODE
          combustionModels
          momentumTransportModels
          fluidThermoMomentumTransportModels
          thermophysicalTransportModels
          rhoReactionThermophysicalTransportModels
          meshTools
          finiteVolume
          radiationModels
          fvOptions
          regionModels
          sampling
          OpenFOAM
          ${EXTRA_LIBS}
)

target_link_libraries(openFOAM-imported ${LIBS})
