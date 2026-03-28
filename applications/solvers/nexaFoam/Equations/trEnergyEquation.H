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
    thermo2T.correctVibEnergy();
}

rhoE.boundaryFieldRef() ==
    rho.boundaryField()*
    (
        hTR.boundaryField() + 0.5*magSqr(U.boundaryField())
    );
