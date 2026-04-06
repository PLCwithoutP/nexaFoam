// First pass using current p
hTR =
    (rhoE + p)/rho
  - 0.5*magSqr(U)
  - eV;

hTR.correctBoundaryConditions();

thermo2T.correct();

// Refresh pressure from updated psi
p.ref() =
    rho()
   /psi();

p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

// Second pass using refreshed p
hTR =
    (rhoE + p)/rho
  - 0.5*magSqr(U)
  - eV;

hTR.correctBoundaryConditions();

thermo2T.correct();

// Keep rhoE BC consistent with conservative total energy
rhoE.boundaryFieldRef() ==
    rho.boundaryField()
   *(
        hTR.boundaryField()
      + eV.boundaryField()
      + 0.5*magSqr(U.boundaryField())
    )
  - p.boundaryField();
  