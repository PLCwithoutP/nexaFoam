qFourierPhi =
    dimensionedScalar
    (
        "zero",
        qFourierPhi.dimensions(),
        0.0
    );
    
if (!inviscid)
{
    volScalarField kappaTRField("kappaTRField", thermo2T.kappaTR());

    if (!use2T)
    {
        qFourierPhi =
        (
            - fvc::interpolate(kappaTRField)
             * mesh.magSf()
             * fvc::snGrad(TTR)
        );
    }
    else
    {
        volScalarField kappaVibField("kappaVibField", thermo2T.kappaVib());

        qFourierPhi =
        (
            - fvc::interpolate(kappaTRField)
             * mesh.magSf()
             * fvc::snGrad(TTR)

            - fvc::interpolate(kappaVibField)
             * mesh.magSf()
             * fvc::snGrad(TVib)
        );
    }
}


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

surfaceScalarField phiJhVib
(
    IOobject
    (
        "phiJhVib",
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

phiJhVib.setOriented(true);

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

if (!inviscid && mixtureCheck && use2T)
{
    forAll(Y_species, speciei)
    {
        if (!vibSpecieActive[speciei])
        {
            continue;
        }

        phiJhVib += phiJYi[speciei]*fvc::interpolate(eVibSpecies[speciei]);
    }
}

solve
(
    fvm::ddt(rhoE)
  + fvc::div(phiEp)
  + fvc::div(qFourierPhi + phiJhTR + phiJhVib)
  - fvc::div(sigmaDotU)
);
