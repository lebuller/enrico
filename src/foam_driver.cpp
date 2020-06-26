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
    args_ = std::make_shared<Foam::argList> (argc, pargv, false, false, true, &comm);
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

    n_fluid_regions_=fluidRegions.size();
    n_solid_regions_=solidRegions.size();
    n_total_regions_=n_fluid_regions_+n_solid_regions_;

    for (int i = 0; i < n_fluid_regions_; i++)
    {
      nelgt_ = nelgt_ + fluidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + fluidRegions[i].nCells();
    }
    for (int i = 0; i < n_solid_regions_; i++)
    {
      nelgt_ = nelgt_ + solidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + solidRegions[i].nCells();
    }

    local_regions_size_.resize(n_total_regions_);
    for (int i = 0; i <n_total_regions_; i++)
    {
      if (i < n_fluid_regions_)
      {
        local_regions_size_.at(i)=fluidRegions[i].nCells();
      }
      else
      {
        local_regions_size_.at(i)=solidRegions[i-n_fluid_regions_].nCells();
      }
    }

    init_displs();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


FoamDriver::initialize() {

}



//std::vector<int> FoamDriver::get_elem(int32_t local_elem)
std::pair<int,int> FoamDriver::get_elem(int32_t local_elem) const
{
  int search_num = 0;
  int temp = 0;

  int region_num;
  int region_elem;

  for (int32_t i = 0; i < n_total_regions_; i++)
  {
    temp = search_num + local_regions_size_.at(i);
    if (temp < local_elem)
    {
      search_num=temp;
    }
    else
    {
      region_num=i;
      region_elem = temp - local_elem;
      break;
    }

  }
  return {region_num, region_elem};
}

std::vector<double> FoamDriver::temperature_local() const
{
  // Each Foam proc finds the temperatures of its local elements
  std::vector<double> local_elem_temperatures(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_elem_temperatures[i] = this->temperature_at(i);
  }

  return local_elem_temperatures;
}

std::vector<int> FoamDriver::fluid_mask_local() const
{
  std::vector<int> local_fluid_mask(nelt_);
  for (int32_t i = 0; i < nelt_; ++i) {
    local_fluid_mask[i] = this->in_fluid_at(i);
  }
  return local_fluid_mask;
}

std::vector<double> FoamDriver::density_local() const
{
  std::vector<double> local_densities(nelt_);

  for (int32_t i = 0; i < nelt_; ++i) {
    if (this->in_fluid_at(i) == 1) {
      auto T = this->temperature_at(i);
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
//   while (pimples.run(runTime)) {
//        #include "readTimeControls.H"
//        #include "readSolidTimeControls.H"
//        #include "compressibleMultiRegionCourantNo.H"
//        #include "solidRegionDiffusionNo.H"
//        #include "setMultiRegionDeltaT.H"

//        runTime++;

//        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
//        while (pimples.loop()) {
//            forAll(fluidRegions, i) {
//!           for (int i = 0; i < fluidRegions.size(); ++i
//                Info<< "\nSolving for fluid region "
//                    << fluidRegions[i].name() << endl;
//                #include "setRegionFluidFields.H"
//                #include "solveFluid.H"
//            }

//            forAll(solidRegions, i) {
//!           for (int i = 0; i < solidRegions.size(); ++i
//                Info<< "\nSolving for solid region "
//                    << solidRegions[i].name() << endl;
//                #include "setRegionSolidFields.H"
//                #include "solveSolid.H"
//            }
//        }

//        runTime.write();

//        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
//            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
//            << nl << endl;
}

Position FoamDriver::centroid_at(int32_t local_elem) const
{
  double x, y, z;
  std::pair<int,int> region;
  region = FoamDriver::get_elem(local_elem);
  if (region.first < n_fluid_regions_) {
    x = fluidRegions[region.first].C()[region.second].component(0);
    y = fluidRegions[region.first].C()[region.second].component(1);
    z = fluidRegions[region.first].C()[region.second].component(2);
  }
  else {
    x = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(0);
    y = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(1);
    z = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(2);
  }
// err_chk(foam_get_local_elem_centroid(local_elem, &x, &y, &z),
// "Could not find centroid of local element " + std::to_string(local_elem));
  return {x, y, z};
}

std::vector<Position> FoamDriver::centroid_local() const
{
  int n_local = this->n_local_elem();
  std::vector<Position> local_element_centroids(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_element_centroids[i] = this->centroid_at(i);
  }
  return local_element_centroids;
}

double FoamDriver::volume_at(int32_t local_elem) const
{
  double volume;
  std::pair<int,int> region;
  region = FoamDriver::get_elem(local_elem);
  if (region.first < n_fluid_regions_) {
    volume = fluidRegions[region.first].V()[region.second];
  }
  else {
    volume = solidRegions[region.first-n_fluid_regions_].V()[region.second];
  }
//  err_chk(foam_get_local_elem_volume(local_elem, &volume),
///          "Could not find volume of local element " + std::to_string(local_elem));
  return volume;
}

std::vector<double> FoamDriver::volume_local() const
{
  int n_local = this->n_local_elem();
  std::vector<double> local_elem_volumes(n_local);
  for (int32_t i = 0; i < n_local; ++i) {
    local_elem_volumes[i] = this->volume_at(i);
  }
  return local_elem_volumes;
}

double FoamDriver::temperature_at(int32_t local_elem) const
{
  double temperature;
  std::pair<int,int> region;
  region = FoamDriver::get_elem(local_elem);
  if (region.first < n_fluid_regions_){
    temperature = thermoFluid[region.first].T()[region.second];
  }
  else {
    temperature = thermos[region.first - n_fluid_regions_].T()[region.second];
  }

  return temperature;
}

int FoamDriver::in_fluid_at(int32_t local_elem) const
{
  //! may can make this more efficient by using n_fluid_regions_ and
  //! total_regions_size
  std::pair<int,int> region;
  region = FoamDriver::get_elem(local_elem);
  if (region.first < n_fluid_regions_){
    return 1;
  }
  else {
    return 0;
  }
}

int FoamDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  Expects(local_elem >= 1 && local_elem <= nelt_);
  std::pair<int,int> region;
  if (region.first < n_fluid_regions_){
    QFluid[region.first].ref()[region.second] = heat;
  }
  else {
    Qsolid[region.first - n_fluid_regions_].ref()[region.second] = heat;
  }
//  return foam_set_heat_source(local_elem, heat);
}

FoamDriver::~FoamDriver()
{
//  delete args_;
}

} // namespace enrico
