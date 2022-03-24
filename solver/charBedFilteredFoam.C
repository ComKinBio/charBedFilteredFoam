/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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
    coalChemistryFoam

Group
    grpLagrangianSolvers

Description
    Transient solver for compressible, turbulent flow, with coal and limestone
    particle clouds, an energy source, and combustion.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulentFluidThermoModel.H"
#include "coalChemistryTurbulenceModel.H"
// #include "basicThermoCloud.H"
#include "coalCloud.H"
#include "psiReactionThermo.H"
#include "CombustionModel.H"
#include "fvOptions.H"
#include "radiationModel.H"
#include "SLGThermo.H"
#include "pimpleControl.H"
#include "simpleControl.H" //added for DBM
#include "pressureControl.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "coarserGrid.H"
#include "gasFilter.H"
#include "sourceFilter.H"
#include "scalarList.H"
#include "scalarIOList.H"

#include "labelList.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

//     #include "addCheckCaseOptions.H"
//     #include "setRootCase.H"
    #include "setRootCaseLists.H" //Tian
    #include "createTime.H"
    #include "createMesh.H"
//     #include "createControl.H" //acture this file is to creat control object according to the head file, here we include two head file, it cause error using this file.
    pimpleControl pimple(mesh);
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "initContinuityErrs.H"

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;
    scalarList timeCounting(5, 0.0);
    scalarList timeFlag(9, 0.0);
    
    Info<< "\ncell volumes\n"<< mesh.V() << endl;
                   
    while (runTime.run())
    {
        if (fixedParticle)
        {
            coalParcels.DEMFlagSet(0);
        }
        else if (useDEMFlag)
        {
        //for DEM 
            if (coalParcels.DEMFlag == 1)
            {
                if (DEMStep>=DEMStepMax)
                {
                    coalParcels.DEMFlagSet(0);
                    DEMStep = 0;
                    DEMNextTime = runTime.value() + DEMInter;
                }
                    
                DEMStep = DEMStep + 1;
            }
            else
            {
                if (runTime.value() >= DEMNextTime)
                {
                    coalParcels.DEMFlagSet(1);
                }
            }
        }
        
        
        
        Info<< "\nDEMFlag\n" <<coalParcels.DEMFlag<< nl << endl;
        
//---
//         Info<< "\nlist\n" <<coalParcels.creatStepList(6,0.060)<< nl << endl;
        
        timeFlag[0] = runTime.elapsedCpuTime();
        
        #include "readTimeControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        pDyn = 0.5*rhoc*magSqr(Uc);

        // Update continuous phase volume fraction field
        alphac = max(1.0 - coalParcels.theta(), alphacMin);
        alphac.correctBoundaryConditions();
        
        timeFlag[1] = runTime.elapsedCpuTime();

        timeFlag[7] = runTime.elapsedCpuTime();
        
        if (useSourceFilter)
        {
            alphacavg = max(1.0 - coalParcels.thetaFiltered(), alphacMin);
            
        }
        
        timeFlag[2] = runTime.elapsedCpuTime();       

        coalParcels.evolve();
        
        timeFlag[3] = runTime.elapsedCpuTime();
        
        if (useGasFilter)
        {
            
            labelListList gasFilterlist = coalParcels.creatStepList(gasFilterStep, gasFilterB, gasFilterRatio);
// Info<< "\nlist\n" <<gasFilterlist<< nl << endl;            
            Ucavg = gasFilterModel.filteredField(Uc, gasFilterlist);
            rhoavg = gasFilterModel.filteredField(rhoc, gasFilterlist);
            Tavg = gasFilterModel.filteredField(T, gasFilterlist);
            muavg = thermo.mu();
            kappaavg = thermo.kappa();
            Cpavg = gasFilterModel.filteredField(thermo.Cp(), gasFilterlist);
            O2avg = gasFilterModel.filteredField(GasYO2, gasFilterlist);
            CO2avg = gasFilterModel.filteredField(GasYCO2, gasFilterlist);
            
        }
        
        timeFlag[8] = runTime.elapsedCpuTime();

        //smooth particle radiation properties
        if (useSourceFilter)
        {
            coalParcels.filerSourceTerms(useSourceFilterStep);
        }
        
        timeFlag[4] = runTime.elapsedCpuTime();
        
        alphacf = fvc::interpolate(alphacavg);
        //phic = linearInterpolate(Uc) & mesh.Sf();
        //rhocPhic = fvc::interpolate(rhoc)*phic;
        alphaRhoPhic = alphacf*rhocPhic;

        fvVectorMatrix averagedSU = coalParcels.SU(Uc);
        
        volScalarField& hee = thermo.he();

        fvScalarMatrix averagedSh = coalParcels.Sh(hee);
        
        timeFlag[5] = runTime.elapsedCpuTime();
        
        #include "rhocEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UcEqn.H"
            #include "YEqn.H"
            #include "EEqn.H"

            //--- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        rhoc = thermo.rho();

        runTime.write();
        
        timeFlag[6] = runTime.elapsedCpuTime();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
            
        //flag 1-7 TFM gas, 3-8 gas filter; 2-3 parcel; 7-2, 8-4 sourcefilter; 4-5 TFM source; 0-1 and 5-6 gas solver
        timeCounting[0]+= timeFlag[6] - timeFlag[0]; //all time
        timeCounting[1]+= timeFlag[1] - timeFlag[0];
        timeCounting[1]+= timeFlag[6] - timeFlag[5]; //gas solver
        timeCounting[2]+= timeFlag[3] - timeFlag[2]; //parcels solver
        timeCounting[3]+= timeFlag[8] - timeFlag[3];
        timeCounting[3]+= timeFlag[7] - timeFlag[1]; //gas average
        timeCounting[4]+= timeFlag[5] - timeFlag[8]; //source average
        
        
        scalarIOList timeStatistic
        (
            IOobject
            (
                    "timeStatistic",
                    runTime.constant(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
            timeCounting
        );
        timeStatistic.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
