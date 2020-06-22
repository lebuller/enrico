
#include "enricoFoamLibrary.H"
#include "fvCFD.H"
#include "fluidThermoMomentumTransportModel.H"
#include "rhoReactionThermophysicalTransportModel.H"
#include "rhoReactionThermo.H"
#include "CombustionModel.H"
#include "fixedGradientFvPatchFields.H"
#include "regionProperties.H"
#include "compressibleCourantNo.H"
#include "solidRegionDiffNo.H"
#include "solidThermo.H"
#include "radiationModel.H"
#include "fvOptions.H"
#include "coordinateSystem.H"
#include "pimpleMultiRegionControl.H"
#include "pressureControl.H"

void foam_init(void *comm)
{

  int argc=2;
  char *argv[1];
  std::string parallel_arg = "-parallel";
  argv[0]=const_cast<char*>( parallel_arg.c_str());
  char** pargv =argv;

  #define NO_CONTROL
  #define CREATE_MESH createMeshesPostProcess.H
// postProcss.H also includes call to setRootCase.H
//  #include "postProcess.H"

// These includes and argList call replace the contents of setRootCaseLists.H
  #include "listOptions.H"
  Foam::argList args(argc, pargv, false, false, true, comm);
  if (!args.checkRootCase())
  {
    Foam::FatalError.exit();
  }
  #include "listOutput.H"

  #include "createTime.H"
  #include "createMeshes.H"
  #include "createFields.H"
  #include "initContinuityErrs.H"
  pimpleMultiRegionControl pimples(fluidRegions, solidRegions);
  #include "createFluidPressureControls.H"
  #include "createTimeControls.H"
  #include "readSolidTimeControls.H"
  #include "compressibleMultiRegionCourantNo.H"
  #include "solidRegionDiffusionNo.H"
  #include "setInitialMultiRegionDeltaT.H"

}
