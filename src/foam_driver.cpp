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
  if (active()) {
    initialize(comm);
    init_displs();
  }
  MPI_Barrier(MPI_COMM_WORLD);
}


void FoamDriver::initialize(MPI_Comm comm) {

    int argc=2;
    char *argv[1];
    std::string parallel_arg = "-parallel";
    argv[0]=const_cast<char*>( parallel_arg.c_str());
    char** pargv =argv;

    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H

    #include "listOptions.H"
    args_ = std::make_shared<Foam::argList> (argc, pargv, false, false, true, &comm);
    Foam::argList &args = *(args_.get());
    if (!args.checkRootCase()) {
      Foam::FatalError.exit();
    }
    #include "listOutput.H"

    #include "createTime.H"

//    #include "createMeshes.H"
    regionProperties rp(runTime);

//      #include "createFluidMeshes.H"
    const wordList fluidNames (rp.found("fluid") ? rp["fluid"] : wordList(0));

    fluidRegions.setSize((fluidNames.size()));

    forAll(fluidNames, i) {
      Info<< "Create fluid mesh for region " << fluidNames[i]
          << " for time = " << runTime.timeName() << nl << endl;

      fluidRegions.set (i, new fvMesh
        (
          IOobject(fluidNames[i], runTime.timeName(), runTime, IOobject::MUST_READ)
        )
      );
    }

//      #include "createsolidMeshes.H"
    const wordList solidNames(rp.found("solid") ? rp["solid"] : wordList(0));

    solidRegions.setSize((solidNames.size()));

    forAll(solidNames, i) {
      Info<< "Create solid mesh for region " << solidNames[i]
          << " for time = " << runTime.timeName() << nl << endl;

      solidRegions.set(i, new fvMesh
        (
          IOobject(solidNames[i], runTime.timeName(), runTime, IOobject::MUST_READ)
        )
      );
    }
//    end #include "createMeshes.H"

//    #include "createFields.H"
//      #include createFluidFields.H
    thermoFluid.setSize(fluidRegions.size());
    rhoFluid.setSize(fluidRegions.size());
    QFluid.setSize(fluidRegions.size());
    UFluid.setSize(fluidRegions.size());
    phiFluid.setSize(fluidRegions.size());
    gFluid.setSize(fluidRegions.size());
    hRefFluid.setSize(fluidRegions.size());
    ghFluid.setSize(fluidRegions.size());
    ghfFluid.setSize(fluidRegions.size());
    turbulenceFluid.setSize(fluidRegions.size());
    thermophysicalTransportFluid.setSize(fluidRegions.size());
    reactionFluid.setSize(fluidRegions.size());
    p_rghFluid.setSize(fluidRegions.size());
    radiation.setSize(fluidRegions.size());
    KFluid.setSize(fluidRegions.size());
    dpdtFluid.setSize(fluidRegions.size());
    fieldsFluid.setSize(fluidRegions.size());
    initialMassFluid.resize(fluidRegions.size());
    MRFfluid.setSize(fluidRegions.size());
    fluidFvOptions.setSize(fluidRegions.size());

// Populate fluid field pointer lists
    forAll(fluidRegions, i) {
      Info<< "*** Reading fluid mesh thermophysical properties for region "
          << fluidRegions[i].name() << nl << endl;
      Info<< "    Adding to thermoFluid\n" << endl;
      thermoFluid.set(i, rhoReactionThermo::New(fluidRegions[i]).ptr());
      Info<< "    Adding to rhoFluid\n" << endl;
      rhoFluid.set (i, new volScalarField
        (
          IOobject("rho", runTime.timeName(), fluidRegions[i], IOobject::NO_READ, IOobject::AUTO_WRITE),
          thermoFluid[i].rho()
        )
      );

      Info<< "    Adding to QFluid\n" << endl;
      QFluid.set(i, new volScalarField
        (
          IOobject("Q", runTime.timeName(), fluidRegions[i], IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE),
          fluidRegions[i]
        )
      );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimLength, 0)
        )
    );

    dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::momentumTransportModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to thermophysicalTransport\n" << endl;
    thermophysicalTransportFluid.set
    (
        i,
        rhoReactionThermophysicalTransportModel::New
        (
            turbulenceFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        CombustionModel<rhoReactionThermo>::New
        (
            thermoFluid[i],
            turbulenceFluid[i]
        )
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Force p_rgh to be consistent with p
    p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    Info<< "    Adding to radiationFluid\n" << endl;
    radiation.set
    (
        i,
        radiationModel::New(thermoFluid[i].T())
    );

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvOptions\n" << endl;
    fluidFvOptions.set
    (
        i,
        new fv::options(fluidRegions[i])
    );

    turbulenceFluid[i].validate();
}

//      #include "createSolidFields.H"
coordinates.setSize(solidRegions.size());
thermos.setSize(solidRegions.size());
radiations.setSize(solidRegions.size());
solidHeatSources.setSize(solidRegions.size());
betavSolid.setSize(solidRegions.size());
Qsolid.setSize(solidRegions.size());
aniAlphas.setSize(solidRegions.size());

forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermos\n" << endl;
    thermos.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to Qsolid\n" << endl;
        Qsolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "Q",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                solidRegions[i]
            )
        );

    Info<< "    Adding to radiations\n" << endl;
    radiations.set(i, radiationModel::New(thermos[i].T()));

    Info<< "    Adding fvOptions\n" << endl;
    solidHeatSources.set
    (
        i,
        new fv::options(solidRegions[i])
    );

    if (!thermos[i].isotropic())
    {
        Info<< "    Adding coordinateSystems\n" << endl;
        coordinates.set
        (
            i,
            coordinateSystem::New(solidRegions[i], thermos[i])
        );

        tmp<volVectorField> tkappaByCp =
            thermos[i].Kappa()/thermos[i].Cp();

        aniAlphas.set
        (
            i,
            new volSymmTensorField
            (
                IOobject
                (
                    "Anialpha",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedSymmTensor
                (
                    "zero",
                    tkappaByCp().dimensions(),
                    Zero
                ),
                zeroGradientFvPatchSymmTensorField::typeName
            )
        );

        aniAlphas[i].primitiveFieldRef() =
            coordinates[i].R().transformVector(tkappaByCp());
        aniAlphas[i].correctBoundaryConditions();

    }

    IOobject betavSolidIO
    (
        "betavSolid",
        runTime.timeName(),
        solidRegions[i],
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    );

    if (betavSolidIO.typeHeaderOk<volScalarField>(true))
    {
        betavSolid.set
        (
            i,
            new volScalarField(betavSolidIO, solidRegions[i])
        );
    }
    else
    {
        betavSolid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "betavSolid",
                    runTime.timeName(),
                    solidRegions[i],
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                solidRegions[i],
                dimensionedScalar(dimless, scalar(1))
            )
        );
    }
}



//    end #include "createFields.H"

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

    for (int i = 0; i < n_fluid_regions_; i++) {
      nelgt_ = nelgt_ + fluidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + fluidRegions[i].nCells();
    }
    for (int i = 0; i < n_solid_regions_; i++) {
      nelgt_ = nelgt_ + solidRegions[i].globalData().nTotalCells();
      nelt_ = nelt_ + solidRegions[i].nCells();
    }

    local_regions_size_.resize(n_total_regions_);
    for (int i = 0; i < n_total_regions_; i++) {
      if (i < n_fluid_regions_) {
        local_regions_size_.at(i)=fluidRegions[i].nCells();
      } else {
        local_regions_size_.at(i)=solidRegions[i-n_fluid_regions_].nCells();
      }
    }

    for (int i=0; i < n_solid_regions_; i++) {
      Pout <<  "volume check " << solidRegions[i].V()[0] << endl;
    }

}



//std::vector<int> FoamDriver::get_elem(int32_t local_elem)
std::pair<int,int> FoamDriver::get_elem(int32_t local_elem) const
{
  int search_num = 0;

  for (int32_t i = 0; i < n_total_regions_; i++) {
    int region_size = local_regions_size_.at(i);
    if (local_elem < search_num + region_size) {
      return {i, local_elem - search_num};
    }
    search_num += region_size;
  }
  throw std::runtime_error{"Value of local element is too high"};
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
  throw std::runtime_error{ "OpenFOAM solve_step not fully implemented" };
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
    x = fluidRegions[region.first].C()[region.second].component(0)*100;
    y = fluidRegions[region.first].C()[region.second].component(1)*100;
    z = fluidRegions[region.first].C()[region.second].component(2)*100;
  } else {
    x = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(0)*100;
    y = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(1)*100;
    z = solidRegions[region.first-n_fluid_regions_].C()[region.second].component(2)*100;
  }
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
    volume = fluidRegions[region.first].V()[region.second]*pow(100,3);
  } else {
    Pout << "proc " << Pstream::myProcNo() << " region "
         << region.first << " element " << region.second << endl;
    Pout << "solid " << region.first-n_fluid_regions_
         << " " << solidRegions[0].C()[0].component(0) << endl;
    Pout << "solid " << region.first-n_fluid_regions_
         << " " << solidRegions[0].V()[0] << endl;
    volume = solidRegions[region.first-n_fluid_regions_].V()[region.second]*pow(100,3);
  }
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
  } else {
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
  } else {
    return 0;
  }
}

int FoamDriver::set_heat_source_at(int32_t local_elem, double heat)
{
  Expects(local_elem >= 1 && local_elem <= nelt_);
  std::pair<int,int> region;
  region = get_elem(local_elem);
  if (region.first < n_fluid_regions_){
    QFluid[region.first].ref()[region.second] = heat;
  } else {
    Qsolid[region.first - n_fluid_regions_].ref()[region.second] = heat;
  }
//  return foam_set_heat_source(local_elem, heat);
}

FoamDriver::~FoamDriver()
{
}

} // namespace enrico
