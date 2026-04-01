qFourierPhi =
    dimensionedScalar
    (
        "zero",
        qFourierPhi.dimensions(),
        0.0
    );

// --- Species diffusive enthalpy flux for Eq. (7)
surfaceScalarField phiJhTR
(
    IOobject
    (
        "phiJhTR",
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

phiJhTR.setOriented(true);

if (!inviscid && mixtureCheck)
{
    forAll(Y_species, speciei)
    {
        volScalarField hsTRi
        (
            IOobject
            (
                "hsTR_" + Y_species[speciei].name(),
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
                hTR.dimensions(),
                0.0
            )
        );

        // Internal field
        forAll(hsTRi, celli)
        {
            hsTRi[celli] = mix.Hs
            (
                speciei,
                p[celli],
                TTR[celli]
            );
        }

        // Boundary field
        forAll(hsTRi.boundaryField(), patchi)
        {
            forAll(hsTRi.boundaryField()[patchi], facei)
            {
                hsTRi.boundaryFieldRef()[patchi][facei] = mix.Hs
                (
                    speciei,
                    p.boundaryField()[patchi][facei],
                    TTR.boundaryField()[patchi][facei]
                );
            }
        }

        phiJhTR += phiJYi[speciei]*fvc::interpolate(hsTRi);
    }
}

solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  + fvc::div(qFourierPhi)
  - fvc::div(sigmaDotU)
);

hTR = rhoE/rho - 0.5*magSqr(U);
hTR.correctBoundaryConditions();
eT.correctBoundaryConditions();
if (use2T)
{
    eR.correctBoundaryConditions();
    eV.correctBoundaryConditions();
}

thermo2T.correct();
thermo2T.correctTEnergy();
if (use2T)
{
    thermo2T.correctREnergy();
}

rhoE.boundaryFieldRef() ==
    rho.boundaryField()*
    (
        hTR.boundaryField() + 0.5*magSqr(U.boundaryField())
    );
