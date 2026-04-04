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


if (!inviscid && use2T)
{
    volScalarField kappaVibField("kappaVibField", thermo2T.kappaVib());

    qVibPhi =
    (
        - fvc::interpolate(kappaVibField)
          * mesh.magSf()
          * fvc::snGrad(TVib)
    );

    qVibPhi.setOriented(true);
}

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

    // ================================================================
    // A) Pure-mixture 2T case
    // ================================================================
    if (!mixtureCheck)
    {
        // No species-based VT/VV/CV source evaluation here
        // unless you later implement a dedicated pure-mixture source model.

        solve
        (
            fvm::ddt(rhoEv)
          + fvc::div(phiEv)
          + fvc::div(qVibPhi)
         ==
            QVzero
        );

        eV = rhoEv/rho;
        eV.correctBoundaryConditions();
    }

    // ================================================================
    // B) Multispecies 2T case
    // ================================================================
    else
    {
        PtrList<volScalarField>& QVT = thermo2T.correctVibSource(TVibSpecies);
        PtrList<volScalarField>& QVV = thermo2T.correctVibVibSource(TVibSpecies);

        PtrList<volScalarField> QCVempty;
        PtrList<volScalarField>* QCVptr = &QCVempty;

        if (useChemistry)
        {
            const PtrList<volScalarField::Internal>& RRchem = chemistryPtr().RR();
            QCVptr = &thermo2T.correctCVSource(RRchem);
        }

        PtrList<volScalarField>& QCV = *QCVptr;

        PtrList<volScalarField>& Y = thermo2T.composition().Y();

        forAll(Y, speciei)
        {
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
                speciei < QVT.size() ? QVT[speciei] : QVzero
            );

            const volScalarField& QVVi =
            (
                speciei < QVV.size() ? QVV[speciei] : QVzero
            );

            const volScalarField& QCVi =
            (
                speciei < QCV.size() ? QCV[speciei] : QVzero
            );

            solve
            (
                fvm::ddt(rhoYEvibSpecies[speciei])
              + fvc::div(phiRhoEvibYi[speciei] + phiJEvib_i)
              + fvc::div(qVibPhi)
             ==
                QVTi + QVVi + QCVi
            );

            rhoYi[speciei].max(SMALL);

            eVibSpecies[speciei] =
                rhoYEvibSpecies[speciei]/rhoYi[speciei];

            eVibSpecies[speciei].correctBoundaryConditions();
        }

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