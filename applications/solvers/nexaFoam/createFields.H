#include "createRDeltaT.H"
#include "vibrationalFieldCheck.H"
autoPtr<ne2TThermo> n2Thermo
(
    ne2TThermo::New(mesh)
);

ne2TThermo& thermo2T = n2Thermo();

volScalarField& hTR = thermo2T.h();
volScalarField& eT = thermo2T.eT();
volScalarField& eR = thermo2T.eR();
autoPtr<volScalarField> eVPtr
(
    new volScalarField 
    (
        IOobject
        (
            "eV",
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
    )
);
volScalarField& eV = eVPtr();
if (vibrationalCheck)
{
    eV = thermo2T.eVib();
    eV.correctBoundaryConditions();
}
/* IOdictionary thermoDict_
(
    IOobject
    (
        "thermophysicalDict",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    )    
);

const scalar& thetaVib_ = thermoDict_.subDict("thermodynamics").get<scalar>("ThetaVibrational"); */

Info<< "Reading field U\n" << endl;
volVectorField U
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

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo2T.rho()
);

volVectorField rhoU
(
    IOobject
    (
        "rhoU",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*U
);

volScalarField rhoE
(
    IOobject
    (
        "rhoE",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho*(hTR + 0.5*magSqr(U))
);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi("phi", fvc::flux(rhoU));

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo2T
    )
);

