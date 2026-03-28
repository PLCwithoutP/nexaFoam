if (mixtureCheck)
{
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    const PtrList<surfaceScalarField>& phiRhoYi = kurganovFlux.phiRhoYi();

    const label nSpec = Y.size();
    const label last = nSpec - 1;

    // --- 1) Prepare Mass Diffusion Flux Term

    forAll(phiJYi, i)
    {
        phiJYi[i] =
            dimensionedScalar
            (
                "zero",
                phiJYi[i].dimensions(),
                0.0
            );
    }

    surfaceScalarField sumPhiJ
    (
        IOobject
        (
            "sumPhiJ",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimMass/dimTime, 0.0)
    );

    sumPhiJ.setOriented(true);

    // Build N-1 diffusive fluxes
    for (label i = 0; i < last; ++i)
    {
        phiJYi[i] = Foam::FickModel::phiJ(rho, Deff, Y[i]);
        sumPhiJ += phiJYi[i];
    }

    // Close the sum exactly
    phiJYi[last] = -sumPhiJ;


   // --- 2) Solve partial densities rhoYi
    for (label i = 0; i < nSpec; ++i)
    {
        solve
        (
            fvm::ddt(rhoYi[i])
        + fvc::div(phiRhoYi[i] + phiJYi[i])
        );
    }

    // --- 3) Build mixture rho from species
    rho = dimensionedScalar("zeroRho", rho.dimensions(), 0.0);

    for (label i = 0; i < nSpec; ++i)
    {
        rho += rhoYi[i];
    }
    rho.max(SMALL);
    rho.correctBoundaryConditions();

    // --- 4) Recover Yi = rhoYi / rho
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;

        // keep bounded
        Y[i].max(0.0);
        Y[i].min(1.0);
        Y[i].correctBoundaryConditions();
    }

    // Renormalize (after bounding)
    volScalarField sumY
    (
        IOobject("sumY", runTime.timeName(), mesh,
                 IOobject::NO_READ, IOobject::NO_WRITE, IOobject::NO_REGISTER),
        mesh,
        dimensionedScalar(dimless, 0.0)
    );
    for (label i = 0; i < nSpec; ++i) sumY += Y[i];

    sumY.max(SMALL);
    for (label i = 0; i < nSpec; ++i) Y[i] /= sumY;

}
