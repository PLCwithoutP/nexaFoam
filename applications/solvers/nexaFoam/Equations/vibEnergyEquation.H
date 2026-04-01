// ************************************************************************* //
// Vibrational energy equation(s)
// ************************************************************************* //

// Always reset placeholder vibrational conductive flux.
// Real model can be filled later.
qVibPhi =
    dimensionedScalar
    (
        "zero",
        qVibPhi.dimensions(),
        0.0
    );

if (use2T)
{
    volScalarField QVzero
    (
        IOobject
        (
            "QVzero",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    PtrList<volScalarField>& QVT = thermo2T.correctVibSource();
    PtrList<volScalarField>& QVV = thermo2T.correctVibVibSource();

    PtrList<volScalarField> QCVempty;
    PtrList<volScalarField>* QCVptr = &QCVempty;

    if (useChemistry)
    {
        const PtrList<volScalarField::Internal>& RRchem = chemistryPtr().RR();
        QCVptr = &thermo2T.correctCVSource(RRchem);
    }

    PtrList<volScalarField>& QCV = *QCVptr;

    volScalarField QVmix
    (
        IOobject
        (
            "QVmix",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh,
        dimensionedScalar("zero", dimEnergy/dimVolume/dimTime, 0.0)
    );

    forAll(QVT, s)
    {
        if (mixtureCheck)
        {
            if (s < vibSpecieActive.size() && vibSpecieActive[s])
            {
                QVmix += QVT[s];
            }
        }
        else
        {
            QVmix += QVT[s];
        }
    }

    forAll(QVV, s)
    {
        if (mixtureCheck)
        {
            if (s < vibSpecieActive.size() && vibSpecieActive[s])
            {
                QVmix += QVV[s];
            }
        }
        else
        {
            QVmix += QVV[s];
        }
    }

    forAll(QCV, s)
    {
        if (mixtureCheck)
        {
            if (s < vibSpecieActive.size() && vibSpecieActive[s])
            {
                QVmix += QCV[s];
            }
        }
        else
        {
            QVmix += QCV[s];
        }
    }

    // ================================================================
    // A) Pure-mixture 2T case
    //    Conservative variable: rhoEv = rho * eV
    //    Convective flux: phiEv
    // ================================================================
    if (!mixtureCheck)
    {
        solve
        (
            fvm::ddt(rhoEv)
            + fvc::div(phiEv)
            + fvc::div(qVibPhi)
            ==
            QVmix
        );

        // Recover specific vibrational energy
        eV = rhoEv/rho;
        eV.correctBoundaryConditions();

        // Vibrational temperature update comes later
        // thermo2T.correctVibEnergy();
    }

    // ================================================================
    // B) Multispecies 2T case
    //    Conservative variables: rhoYEvibSpecies[i] = rho*Yi*eVib_i
    //    Convective fluxes: phiRhoEvibYi[i]
    //    Diffusive fluxes:  phiJYi[i] * interpolate(eVibSpecies[i])
    // ================================================================
    else
    {
        PtrList<volScalarField>& Y = thermo2T.composition().Y();
        const label nSpec = Y.size();

        forAll(Y, speciei)
        {
            // Skip atomic species completely
            if (!vibSpecieActive[speciei])
            {
                rhoYEvibSpecies[speciei] =
                    dimensionedScalar
                    (
                        "zero",
                        rhoYEvibSpecies[speciei].dimensions(),
                        0.0
                    );

                eVibSpecies[speciei] =
                    dimensionedScalar
                    (
                        "zero",
                        eVibSpecies[speciei].dimensions(),
                        0.0
                    );

                TVibSpecies[speciei] =
                    dimensionedScalar
                    (
                        "zero",
                        TVibSpecies[speciei].dimensions(),
                        0.0
                    );

                rhoYEvibSpecies[speciei].correctBoundaryConditions();
                eVibSpecies[speciei].correctBoundaryConditions();
                TVibSpecies[speciei].correctBoundaryConditions();
                continue;
            }

            // Species diffusive vibrational-energy flux:
            // J_i * eVib_i
            surfaceScalarField phiJEvib_i
            (
                IOobject
                (
                    "phiJEvib_" + Y[speciei].name(),
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    IOobject::NO_REGISTER
                ),
                phiJYi[speciei]*fvc::interpolate(eVibSpecies[speciei])
            );

            phiJEvib_i.setOriented(true);

            const volScalarField& QVTi =
            (
                speciei < QVT.size()
                ? QVT[speciei]
                : QVzero
            );

            const volScalarField& QVVi =
            (
                speciei < QVV.size()
                ? QVV[speciei]
                : QVzero
            );

            const volScalarField& QCVi =
            (
                speciei < QCV.size()
                ? QCV[speciei]
                : QVzero
            );

            solve
            (
                fvm::ddt(rhoYEvibSpecies[speciei])
                + fvc::div(phiRhoEvibYi[speciei] + phiJEvib_i)
                + fvc::div(qVibPhi)
                ==
                QVTi + QVVi + QCVi
            );

            // Recover specific species vibrational energy
            rhoYi[speciei].max(SMALL);

            eVibSpecies[speciei] =
                rhoYEvibSpecies[speciei]/rhoYi[speciei];

            eVibSpecies[speciei].correctBoundaryConditions();

            // Species vibrational temperature update comes later
            // TVibSpecies[speciei] = Newton-Raphson inversion from eVibSpecies[speciei]
        }

        // Reconstruct mixture vibrational energy for post-processing / compatibility
        eV = dimensionedScalar("zero", eV.dimensions(), 0.0);

        forAll(Y, speciei)
        {
            if (!vibSpecieActive[speciei])
            {
                continue;
            }

            eV += rhoYEvibSpecies[speciei]/rho;
        }

        eV.correctBoundaryConditions();
        rhoEv = rho*eV;
        rhoEv.correctBoundaryConditions();
    }
}