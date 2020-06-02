/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    chtMultiRegionFoam

Description
    Solver for steady or transient fluid flow and solid heat conduction, with
    conjugate heat transfer between regions, buoyancy effects, turbulence,
    reactions and radiation modelling.

\*---------------------------------------------------------------------------*/

#include <fstream>
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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #define CREATE_MESH createMeshesPostProcess.H
    #include "postProcess.H"

    #include "setRootCaseLists.H"
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


    int n_local_elements=0;
    int cellI;
//    int procID;
    int np=Pstream::nProcs();
    int n=np-1;

    if(Pstream::myProcNo()==n){
        forAll(fluidRegions, i)
        {
            n_local_elements=n_local_elements+fluidRegions[i].nCells() ;
            Pout<< "\n fluid " << i << " gives " << n_local_elements << " adding "
                << fluidRegions[i].nCells() << " on processor " <<Pstream::myProcNo() <<endl;
        }
        forAll(solidRegions, i)
        {
            n_local_elements=n_local_elements+solidRegions[i].nCells() ;
            Pout<< "\n solid " << i << " gives " << n_local_elements << " adding "
                << solidRegions[i].nCells() << " on processor " << Pstream::myProcNo() <<endl;
        }
        Pout<< "\n The number of local elements is on " << Pstream::myProcNo() << " is "
            << n_local_elements << endl;
    }


    if(Pstream::master())
    {
        std::ofstream volfile;
        std::ofstream centfile;

        string volfile_name;
        string centfile_name;

        forAll(fluidRegions, i)
        {
            volfile_name="cellVolumes"+fluidNames[i]+std::to_string(Pstream::myProcNo())+".txt";
            centfile_name="cellCenters"+fluidNames[i]+std::to_string(Pstream::myProcNo())+".txt";

            volfile.open(volfile_name);
            centfile.open(centfile_name);

            for (cellI=0; cellI<fluidRegions[i].nCells(); cellI++)
            {
                volfile  << "Cell volume " << cellI << " is " << fluidRegions[i].V()[cellI] << "\n";
                centfile << "Cell center " << cellI << " is "
                         << fluidRegions[i].C()[cellI].component(0) << ", "
                         << fluidRegions[i].C()[cellI].component(1) << ", "
                         << fluidRegions[i].C()[cellI].component(2) << "\n";
            }

            volfile.close();
            centfile.close();
        }
        forAll(solidRegions, i)
        {
            volfile_name="cellVolumes"+solidNames[i]+std::to_string(Pstream::myProcNo())+".txt";
            centfile_name="cellCenters"+solidNames[i]+std::to_string(Pstream::myProcNo())+".txt";

            volfile.open(volfile_name);
            centfile.open(centfile_name);
            for (cellI=0; cellI<solidRegions[i].nCells(); cellI++){
                volfile  << "Cell volume " << cellI << " is " << solidRegions[i].V()[cellI] << "\n";
                centfile << "Cell center " << cellI << " is "
                         << solidRegions[i].C()[cellI].component(0) << ", "
                         << solidRegions[i].C()[cellI].component(1) << ", "
                         << solidRegions[i].C()[cellI].component(2) << "\n";
            }

            volfile.close();
            centfile.close();
        }
    }

    while (pimples.run(runTime))
    {
        #include "readTimeControls.H"
        #include "readSolidTimeControls.H"

        #include "compressibleMultiRegionCourantNo.H"
        #include "solidRegionDiffusionNo.H"
        #include "setMultiRegionDeltaT.H"

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- PIMPLE loop
        while (pimples.loop())
        {
            forAll(fluidRegions, i)
            {
                Info<< "\nSolving for fluid region "
                    << fluidRegions[i].name() << endl;
                #include "setRegionFluidFields.H"
                #include "solveFluid.H"
            }

            forAll(solidRegions, i)
            {
                Info<< "\nSolving for solid region "
                    << solidRegions[i].name() << endl;
                #include "setRegionSolidFields.H"
                #include "solveSolid.H"
            }
        }

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }


    if(Pstream::master())
    {
        std::ofstream tfile;
        std::ofstream rhofile;

        string tfile_name;
        string rhofile_name;

        forAll(fluidRegions, i)
        {
            tfile_name="cellTemperature"+fluidNames[i]+std::to_string(Pstream::myProcNo())+".txt";
            rhofile_name="cellDensity"+fluidNames[i]+std::to_string(Pstream::myProcNo())+".txt";

            tfile.open(tfile_name);
            rhofile.open(rhofile_name);

            for (cellI=0; cellI<fluidRegions[i].nCells(); cellI++)
            {
                tfile   << "Cell temperature " << cellI << " is "
                        << thermoFluid[i].T()[cellI] << "\n";
                rhofile << "Cell density "     << cellI << " is "
                        << thermoFluid[i].rho()[cellI]        << "\n";
            }

            tfile.close();
            rhofile.close();
        }
        forAll(solidRegions, i)
        {
            tfile_name="cellTemperature"+solidNames[i]+std::to_string(Pstream::myProcNo())+".txt";
            rhofile_name="cellTemperature"+solidNames[i]+std::to_string(Pstream::myProcNo())+".txt";

            tfile.open(tfile_name);
            rhofile.open(rhofile_name);

            for (cellI=0; cellI<solidRegions[i].nCells(); cellI++)
            {
                tfile   << "Cell temperature " << cellI << " is "
                        << thermos[i].T()[cellI] << "\n";
                rhofile   << "Cell density " << cellI << " is "
                        << thermos[i].rho()[cellI] << "\n";
            }

            tfile.close();
            rhofile.close();
        }
    }


    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
