#include "createRDeltaT.H"

#include "readGravitationalAcceleration.H"

Info<< "Reading thermophysical properties\n" << endl;
autoPtr<psiReactionThermo> pThermo(psiReactionThermo::New(mesh));
psiReactionThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

SLGThermo slgThermo(mesh, thermo);

basicSpecieMixture& composition = thermo.composition();
PtrList<volScalarField>& Y = composition.Y();

label O2Id = slgThermo.carrierId("O2");
label CO2Id = slgThermo.carrierId("CO2");
label COId = slgThermo.carrierId("CO");

volScalarField& GasYO2 = composition.Y(O2Id);
volScalarField& GasYCO2 = composition.Y(CO2Id);

volScalarField& p = thermo.p();
const volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField O2avg
(
    IOobject
    (
        "O2avg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    GasYO2
);

volScalarField H2Oavg
(
    IOobject
    (
        "H2Oavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    GasYO2
);

volScalarField CO2avg
(
    IOobject
    (
        "CO2avg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    GasYCO2
);

volScalarField H2avg
(
    IOobject
    (
        "H2avg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    GasYO2
);

const word inertSpecie(thermo.lookup("inertSpecie"));
if (!composition.species().found(inertSpecie))
{
    FatalIOErrorIn(args.executable().c_str(), thermo)
        << "Inert specie " << inertSpecie << " not found in available species "
        << composition.species()
        << exit(FatalIOError);
}

Info<< "inertSpecie" <<inertSpecie<< nl << endl;

volScalarField Tavg
(
    IOobject
    (
        "Tavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    T
);

volScalarField kappaavg
(
    IOobject
    (
        "kappaavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.kappa()
);
    
volScalarField Cpavg
(
    IOobject
    (
        "cpavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.Cp()
);

volScalarField muavg
(
    IOobject
    (
        "muavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.mu()
);


multivariateSurfaceInterpolationScheme<scalar>::fieldTable fields;

forAll(Y, i)
{
    fields.add(Y[i]);
}
fields.add(thermo.he());

volScalarField rhoc
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField rhoavg
(
    IOobject
    (
        "rhoavg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

// dynamic pressure field - used externally (optional)
volScalarField pDyn
(
    IOobject
    (
        "pDyn",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("zero",dimPressure, 0.0)
);

//change --
Info<< "\nReading field Uc\n" << endl;
volVectorField Uc
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volVectorField Ucavg
(
    IOobject
    (
        "Uavg",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Uc
);
//--
//change --
//#include "compressibleCreatePhi.H"
Info<< "Reading/calculating face flux field phic\n" << endl;

surfaceScalarField phic
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(Uc) & mesh.Sf()
);

label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, pimple.dict(), pRefCell, pRefValue);

mesh.setFluxRequired(p.name());

volScalarField Qdot
(
    IOobject
    (
        "Qdot",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Qdot",dimEnergy/dimVolume/dimTime, 0.0)
);


volScalarField alphac
(
    IOobject
    (
        "alphac",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("alphac",dimless, 0.0)
);

volScalarField alphacavg
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphac
);

// #include "createDpdt.H"
Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);//Tian

//change --
//#include "createK.H"
Info<< "Creating field kinetic energy Kc\n" << endl;
volScalarField Kc("K", 0.5*magSqr(Uc));

if (Uc.nOldTimes())
{
    volVectorField* Ucold = &Uc.oldTime();
    volScalarField* Kcold = &Kc.oldTime();
    *Kcold == 0.5*magSqr(*Ucold);

    while (Ucold->nOldTimes())
    {
        Ucold = &Ucold->oldTime();
        Kcold = &Kcold->oldTime();
        *Kcold == 0.5*magSqr(*Ucold);
    }
}

#include "createMRF.H"

#include "filter.H"

#include "createClouds.H"

volScalarField parcelNumber
(
    IOobject
    (
        "parcelNumber",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    coalParcels.np()
);

volScalarField parcelVolumeField
(
    IOobject
    (
        "parcelVolumeField",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    coalParcels.parcelVolumeField()
);

scalar alphacMin
(
    readScalar
    (
        coalParcels.particleProperties().subDict("constantProperties")
       .lookup("alphaMin")
    )
);

// Update alphac from the particle locations
alphac = max(1.0 - coalParcels.theta(), alphacMin);
alphac.correctBoundaryConditions();

surfaceScalarField alphacf("alphaf", fvc::interpolate(alphacavg));
// Info<<"alphacf: "<<alphacf<<endl;
surfaceScalarField rhocPhic
(
    IOobject
    (
        "rhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    fvc::interpolate(rhoc)*phic
);

surfaceScalarField alphaRhoPhic
(
    IOobject
    (
        "alphaRhoPhi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    alphacf*rhocPhic
);


Info<< "Creating turbulence model\n" << endl;
autoPtr<coalChemistryTurbulenceModel> turbulence
(
    coalChemistryTurbulenceModel::New
    (
        alphacavg,
        rhoc,
        Uc,
        alphaRhoPhic,
        phic,
        thermo
    )
);

Info<< "Creating combustion model\n" << endl;
autoPtr<CombustionModel<psiReactionThermo>> combustion
(
    CombustionModel<psiReactionThermo>::New(thermo, turbulence())
);


//--
#include "createRadiationModel.H"

#include "createFvOptions.H"

volScalarField Gavg
(
    IOobject
    (
        "Gavg",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Gavg",dimMass/pow3(dimTime), 0.0)
);




