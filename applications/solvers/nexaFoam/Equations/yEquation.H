if (mixtureCheck)
{
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    const PtrList<surfaceScalarField>& phiRhoYi = kurganovFlux.phiRhoYi();

    const label nSpec = Y.size();

    // --- 1) Solve partial densities rhoYi
    for (label i = 0; i < nSpec; ++i)
    {
        solve
        (
            fvm::ddt(rhoYi[i]) + fvc::div(phiRhoYi[i])
        );
    }

    // --- 2) Build mixture rho from species
    rho = dimensionedScalar("zeroRho", rho.dimensions(), 0.0);

    for (label i = 0; i < nSpec; ++i)
    {
        rho += rhoYi[i];
    }
    rho.max(SMALL);
    rho.correctBoundaryConditions();

    // --- 3) Recover Yi = rhoYi / rho
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;

        // keep bounded
        Y[i].max(0.0);
        Y[i].min(1.0);
        Y[i].correctBoundaryConditions();
    }

    // Optional renormalize (after bounding)
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

// Sanity check
if (mixtureCheck)
{
    volScalarField drho(rhoCont - rho);
    const scalar maxErr = gMax(mag(drho)().primitiveField());
    Info<< "Sanity check: max(|rhoCont - rhoFromSpecies|) = "
        << maxErr << nl << endl;
}
