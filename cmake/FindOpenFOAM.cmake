if(NOT DEFINED ENV{WM_PROJECT})
  message(WARNING "Environment does not contain the necessary OpenFOAM variables.")
  set(OpenFOAM_FOUND False)
else()

  set(OPENFOAM_LINK_DIRS $ENV{FOAM_LIBBIN})

  # removed -m64 flag from this list
  set(OPENFOAM_DEFINITIONS -Dlinux64 -DWM_ARCH_OPTION=64 -DWM_DP -DWM_LABEL_SIZE=32 -DNoRepository)
  set(OPENFOAM_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Xlinker --no-as-needed -Xlinker --add-needed")

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
  $ENV{FOAM_SRC}/meshTools/lnInclude/
  $ENV{FOAM_SRC}/sampling/lnInclude/
  $ENV{FOAM_SRC}/finiteVolume/lnInclude/
  $ENV{FOAM_SRC}/finiteVolume/cfdTools/
  $ENV{FOAM_SRC}/thermophysicalModels/basic/lnInclude/
  $ENV{FOAM_SRC}/thermophysicalModels/specie/lnInclude/
  $ENV{FOAM_SRC}/thermophysicalModels/reactionThermo/lnInclude/
  $ENV{FOAM_SRC}/thermophysicalModels/solidThermo/lnInclude/
  $ENV{FOAM_SRC}/thermophysicalModels/chemistryModel/lnInclude/
  $ENV{FOAM_SRC}/ODE/lnInclude/
  $ENV{FOAM_SRC}/combustionModels/lnInclude/
  $ENV{FOAM_SRC}/MomentumTransportModels/momentumTransportModels/lnInclude/
  $ENV{FOAM_SRC}/MomentumTransportModels/compressible/lnInclude/
  $ENV{FOAM_SRC}/ThermophysicalTransportModels/lnInclude/
  $ENV{FOAM_SRC}/ThermophysicalTransportModels/rhoReactionThermo/lnInclude/
  $ENV{FOAM_SRC}/radiationModels/lnInclude/
  $ENV{FOAM_SRC}/regionModels/regionModel/lnInclude/
  $ENV{FOAM_SRC}/OpenFOAM/lnInclude/
  $ENV{FOAM_SRC}/OSspecific/POSIX/lnInclude/
  )

  set(OPENFOAM_EXTRA_LIBS dl m)
  set(OPENFOAM_LIBRARY_NAMES
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

  # get absolute paths of the necessary librares
  foreach(FOAM_LIB ${OPENFOAM_LIBRARY_NAMES})
    find_library("${FOAM_LIB}_LIB" ${FOAM_LIB} HINTS ${OPENFOAM_LINK_DIRS})
    list(APPEND OPENFOAM_LIBRARIES ${${FOAM_LIB}_LIB})
  endforeach(FOAM_LIB ${OPENFOAM_LIBRARY_NAMES})

#  add_executable(chtMultiRegionFoam ${CHT_DIR}/chtMultiRegionFoam.C ${OPENFOAM_SOURCES})
#  target_compile_definitions(chtMultiRegionFoam PUBLIC ${OPENFOAM_DEFINITIONS})
#  set_target_properties(chtMultiRegionFoam PROPERTIES LINK_FLAGS ${OPENFOAM_LINKER_FLAGS})
#  target_link_libraries(chtMultiRegionFoam PUBLIC ${OPENFOAM_LIBRARIES})
#  target_include_directories(chtMultiRegionFoam PUBLIC ${OPENFOAM_INCLUDE_DIRS})

  ### FOR FUTURE USE, IMPORTED TARGET FOR OPENFOAM ###
  add_library(openFOAM-imported ${CHT_DIR}/enricoFoamLibrary.C ${OPENFOAM_SOURCES})
  target_compile_definitions(openFOAM-imported PUBLIC ${OPENFOAM_DEFINITIONS})
  set_target_properties(openFOAM-imported PROPERTIES LINK_FLAGS ${OPENFOAM_LINKER_FLAGS})
  target_link_libraries(openFOAM-imported PUBLIC ${OPENFOAM_LIBRARIES})
  target_include_directories(openFOAM-imported PUBLIC ${OPENFOAM_INCLUDE_DIRS})

  set(OpenFOAM_FOUND True)
endif()
