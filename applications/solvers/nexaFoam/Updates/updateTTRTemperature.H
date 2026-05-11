//Info<< "Entering updateTTRTemperature" << nl << endl;

thermo2T.eVib() = rhoEv / rho;

eTR =
    rhoE / rho
  - 0.5*magSqr(U)
  - eV;

eTR.correctBoundaryConditions();  // required for mixed BC
TTR.correctBoundaryConditions();  // required for mixed BC
thermo2T.correct();               

// Keep rhoE BC consistent with conservative total energy
rhoE.boundaryFieldRef() ==
    rho.boundaryField()
   *(
        eTR.boundaryField()
      + eV.boundaryField()
      + 0.5*magSqr(U.boundaryField())
    );

// Refresh pressure from updated psi
p.ref() = rho() / psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
