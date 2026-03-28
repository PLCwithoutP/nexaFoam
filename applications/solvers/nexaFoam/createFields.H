#include "createRDeltaT.H"
#include "thermochemistryCheck.H"
#include "createMixture.H"

volScalarField& p = thermo2T.p();
const volScalarField& TTR = thermo2T.TTR();
const volScalarField& TVib = thermo2T.TVib();
const volScalarField& psi = thermo2T.psi();

bool inviscid(true);
if (max(thermo2T.mu().cref().primitiveField()) > 0.0)
{
    inviscid = false;
}

PtrList<volScalarField>& Y_species = mix.Y();

volScalarField& hTR = thermo2T.h();
volScalarField& eT = thermo2T.eT();

autoPtr<volScalarField> eRPtr
(
    new volScalarField 
    (
        IOobject
        (
            "eR",
            runTime.timeName(), 
            mesh, 
            IOobject::NO_READ, 
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimMass, 0.0)
    )
);
volScalarField& eR = eRPtr();

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

if (use2T)
{
    eR = thermo2T.eR();
    eR.correctBoundaryConditions();

    eV = thermo2T.eVib();
    eV.correctBoundaryConditions();
}
else
{
    eR = dimensionedScalar("zero", dimEnergy/dimMass, 0.0);
    eV = dimensionedScalar("zero", dimEnergy/dimMass, 0.0);

    eR.correctBoundaryConditions();
    eV.correctBoundaryConditions();
}

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

Info<< "Turbulence model: Laminar\n" << endl;


// Mixture rho for sanity check purposes
volScalarField rhoCont
(
    IOobject("rhoCont", runTime.timeName(), mesh,
             IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
    rho
);

// Mixture fluxes
PtrList<volScalarField> rhoYi;

if (mixtureCheck)
{
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    rhoYi.setSize(Y.size());

    forAll(Y, i)
    {
        const word name("rhoY_" + Y[i].name());

        rhoYi.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    name,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                rho*Y[i]
            )
        );
    }
}

// Effective diffusion coefficient field [m^2/s]
volScalarField Deff
(
    IOobject
    (
        "Deff",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE     
    ),
    mesh,
    dimensionedScalar
    (
        "zero", 
        dimArea/dimTime, 
        0.0
    )
);

// --- Fourier conductive heat flux for Eq. (7): Translational-Rotational Energy Equation
surfaceScalarField qFourierPhi
(
    IOobject
    (
        "qFourierPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimEnergy/dimTime,
        0.0
    )
);

PtrList<surfaceScalarField> phiJYi;

if (mixtureCheck)
{
    const PtrList<volScalarField>& Y = thermo2T.composition().Y();

    phiJYi.setSize(Y.size());

    forAll(Y, i)
    {
        phiJYi.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phiJYi_" + Y[i].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimMass/dimTime,
                    0.0
                )
            )
        );

        phiJYi[i].setOriented(true);
    }
}
