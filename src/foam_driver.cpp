#include "enrico/foam_driver.h"
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


#include "enrico/error.h"
#include "gsl/gsl"
#include "iapws/iapws.h"

#include "xtensor/xadapt.hpp"

#include <climits>
#include <string>
#include <unistd.h>
#include <vector>

namespace enrico {

FoamDriver::FoamDriver(MPI_Comm comm, pugi::xml_node node)
  : HeatFluidsDriver(comm,node)
{
  //! Need to determine what this portion of the nekDriver does and how it needs
  //! to be applied for an OpenFOAM driver and setting up MPI
  if (active()) {

    int argc=2;
    char *argv[1];
    std::string parallel_arg = "-parallel";
    argv[0]=const_cast<char*>( parallel_arg.c_str());
    char** pargv =argv;

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H

// These includes and argList call replace the contents of setRootCaseLists.H
    #include "listOptions.H"
//    args_ = new Foam::argList(argc,pargv, false, false, true, (void*)&comm);
    args_ = std::make_shared<Foam::argList> (argc, pargv, false, false, true, (void*)&comm);
    Foam::argList &args = *(args_.get());
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

    nelt_ = 0;

    nelgt_ = 0;

    for (int i = 0; i < fluidRegions.size(); i++)
    {
      nelgt_ = nelgt_ + fluidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + fluidRegions[i].nCells();
    }
    for (int i = 0; i < solidRegions.size(); i++)
    {
      nelgt_ = nelgt_ + solidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + solidRegions[i].nCells();
    }
    init_displs();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

std::vector<double> FoamDriver::temperature_local() const
{
  // Each Foam proc finds the temperatures of its local elements
  std::vector<double> local_elem_temperatures(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_elem_temperatures[i] = this->temperature_at(i + 1);
  }

  return local_elem_temperatures;
}

std::vector<int> FoamDriver::fluid_mask_local() const
{
  std::vector<int> local_fluid_mask(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_fluid_mask[i] = this->in_fluid_at(i + 1);
  }
  return local_fluid_mask;
}

std::vector<double> FoamDriver::density_local() const
{
  std::vector<double> local_densities(nelt_);

  for (int32_t i = 0; i < nelt_; ++i) {
    if (this->in_fluid_at(i + 1) == 1) {
      auto T = this->temperature_at(i + 1);
      // nu1 returns specific volume in [m^3/kg]
      local_densities[i] = 1.0e-3 / iapws::nu1(pressure_bc_, T);
    } else {
      local_densities[i] = 0.0;
    }
  }

  return local_densities;
}

void FoamDriver::solve_step()
{
  //! add calls to the foam module for running x timesteps
}

Position FoamDriver::centroid_at(int32_t local_elem) const
{
  double x, y, z;
//  err_chk(foam_get_local_elem_centroid(local_elem, &x, &y, &z),
//          "Could not find centroid of local element " + std::to_string(local_elem));
  return {x, y, z};
}

std::vector<Position> FoamDriver::centroid_local() const
{
  int n_local = this->n_local_elem();
  std::vector<Position> local_element_centroids(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_element_centroids[i] = this->centroid_at(i + 1);
  }
  return local_element_centroids;
}

double FoamDriver::volume_at(int32_t local_elem) const
{
  double volume;
//  err_chk(foam_get_local_elem_volume(local_elem, &volume),
///          "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

std::vector<double> FoamDriver::volume_local() const
{
  int n_local = this->n_local_elem();
  std::vector<double> local_elem_volumes(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_elem_volumes[i] = this->volume_at(i + 1);
  }
  return local_elem_volumes;
}

double FoamDriver::temperature_at(int32_t local_elem) const
{
  double temperature;
//  err_chk(foam_get_local_elem_temperature(local_elem, &temperature),
//          "Could not find temperature of local element " + std::to_string(local_elem));
  return temperature;
}

int FoamDriver::in_fluid_at(int32_t local_elem) const
{
  //! This can potentially be (May need to be) wrapped into the init portion that gets
  //! number of elements locally and globally (which also need to be different than Nek)
//  return foam_local_elem_is_in_fluid(local_elem);
}

int FoamDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  Expects(local_elem >= 1 && local_elem <= nelt_);
//  return foam_set_heat_source(local_elem, heat);
}

FoamDriver::~FoamDriver()
{
//  delete args_;
}

} // namespace enrico
