if (useChemistry)
{
    const scalar dtChem = runTime.deltaTValue();

    // 1) Solve the local chemistry ODE system.
    //    This internally uses the chemistry solver selected in chemistryProperties (currently neEulerImplicit).
    chemistryPtr->solve(dtChem);

    // 2) Species source-term split update: rhoYi^{n+1} = rhoYi^* + dt * RR_i
    PtrList<volScalarField>& Y = thermo2T.composition().Y();
    const label nSpec = Y.size();

    for (label i = 0; i < nSpec; ++i)
    {
        const volScalarField::Internal& RRi = chemistryPtr->RR(i);

        forAll(rhoYi[i], celli)
        {
            rhoYi[i][celli] += dtChem*RRi[celli];

            if (rhoYi[i][celli] < SMALL)
            {
                rhoYi[i][celli] = SMALL;
            }
        }

        rhoYi[i].correctBoundaryConditions();
    }

    // 3) Rebuild mixture density from partial densities
    rho = dimensionedScalar("zeroRho", rho.dimensions(), 0.0);

    for (label i = 0; i < nSpec; ++i)
    {
        rho += rhoYi[i];
    }

    rho.max(SMALL);
    rho.correctBoundaryConditions();

    // 4) Recover mass fractions Yi = rhoYi / rho
    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] = rhoYi[i]/rho;
        Y[i].max(0.0);
        Y[i].min(1.0);
        Y[i].correctBoundaryConditions();
    }

    // 5) Renormalize Y
    volScalarField sumY
    (
        IOobject
        (
            "sumYChemistry",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar(dimless, 0.0)
    );

    for (label i = 0; i < nSpec; ++i)
    {
        sumY += Y[i];
    }

    sumY.max(SMALL);

    for (label i = 0; i < nSpec; ++i)
    {
        Y[i] /= sumY;
        Y[i].correctBoundaryConditions();
    }

    // 6) Apply chemistry heat source to total energy density
    tmp<volScalarField> tQdot = chemistryPtr->Qdot();
    const volScalarField& Qdot = tQdot();

    forAll(rhoE, celli)
    {
        rhoE[celli] += dtChem*Qdot[celli];
    }

    rhoE.correctBoundaryConditions();

    // 7) Rebuild sensible enthalpy/internal-energy variable used by thermo
    hTR =
        (rhoE + p)/rho
        - 0.5*magSqr(U)
        - eV;
    hTR.correctBoundaryConditions();

    // Keep momentum consistent with unchanged velocity
    rhoU = rho*U;
    rhoU.correctBoundaryConditions();

    // 8) Re-correct thermo from the chemistry-updated state
    thermo2T.correct();

    if (use2T)
    {
        thermo2T.correctVibEnergy();
    }

    // 9) Rebuild pressure from updated rho and psi
    p.ref() = rho()/psi();
    p.correctBoundaryConditions();

    // Optional consistency for mixture-check path
    if (mixtureCheck)
    {
        rhoCont = rho;
    }
}
