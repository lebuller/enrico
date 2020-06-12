

if(NOT DEFINED ENV{WM_PROJECT})
  message(WARNING "Environment does not contain the necessary OpenFOAM variables.")
  set(OpenFOAM_FOUND False)
else()
  set(OPENFOAM_LINK_DIRS $ENV{FOAM_LIBBIN})

  # removed -m64 flag from this list
  set(OPENFOAM_DEFINITIONS -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository)

  set(CHT_DIR ${CMAKE_SOURCE_DIR}/cht_interface)

  set(OPENFOAM_SOURCES
  ${CHT_DIR}/solid/solidRegionDiffNo.C
  ${CHT_DIR}/fluid/compressibleCourantNo.C)

  set(OPENFOAM_INCLUDE_DIRS
  ${CHT_DIR}/
  ${CHT_DIR}/fluid/
  ${CHT_DIR}/solid/
  ${CHT_DIR}/porousFluid/
  ${CHT_DIR}/porousSolid/
  ${CHT_DIR}/include/
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

  set(OPENFOAM_EXTRA_LIBS dl m)
  set(OPENFOAM_LIBS
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
  ${OPENFOAM_EXTRA_LIBS}
  )

  add_executable(chtMultiRegionFoam ${CHT_DIR}/chtMultiRegionFoam.C ${OPENFOAM_SOURCES})
  target_compile_definitions(chtMultiRegionFoam PUBLIC ${OPENFOAM_DEFINITIONS})
  target_link_libraries(chtMultiRegionFoam PUBLIC ${OPENFOAM_LIBS})
  # necessary for CMake < 3.12
  # target_link_directories(chtMultiRegionFoam ${OPENFOAM_LINK_DIRS})
  set_target_properties(chtMultiRegionFoam PROPERTIES
  INTERFACE_LINK_DIRECTORIES ${OPENFOAM_LINK_DIRS})
  target_include_directories(chtMultiRegionFoam PUBLIC ${OPENFOAM_INCLUDE_DIRS})

  # add_library(openFOAM-imported INTERFACE IMPORTED)
  # target_link_libraries(openFOAM-imported ${OPENFOAM_LIBS})
  # # necessary for CMake < 3.12
  # # target_link_directories(openFOAM-imported ${OPENFOAM_LINK_DIRS})
  # set_target_properties(openFOAM-imported PROPERTIES
  # INTERFACE_LINK_DIRECTORIES ${OPENFOAM_LINK_DIRS})
  # target_include_directories(openFOAM-imported ${OPENFOAM_INCLUDE_DIRS})

  set(OpenFOAM_FOUND True)
endif()
