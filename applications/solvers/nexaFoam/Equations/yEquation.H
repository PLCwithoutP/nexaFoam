/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Specie continuity equation implementation (no diffusion yet).
    - Active only for multiComponent2TMixture (mixtureCheck == true)
    - Uses Kurganov convective fluxes: phiRhoYi from kurganovFlux
    - Always keeps mixture continuity (rho) equation separate (do not remove).

    Notes:
    - Solves N-1 species and closes with last species: Y_last = 1 - sum(Y_i)
    - Adds basic bounding to keep Yi in [0,1]
\*---------------------------------------------------------------------------*/

if (mixtureCheck)
{
    // Mass-fraction fields from thermo mixture
    PtrList<volScalarField>& Y = thermo2T.composition().Y();

    // Convective species mass fluxes (kg/s) computed inside Kurganov
    const PtrList<surfaceScalarField>& phiRhoYi = kurganovFlux.phiRhoYi();

    const label nSpecie = Y.size();

    if (nSpecie < 2)
    {
        // Nothing to solve
        Info << nl << "yEquation.H: mixtureCheck=true but Y.size()<2. Skipping species solve." << nl << endl;
    }
    else
    {

        // Solve N-1 species to enforce sum(Y)=1
        const label last = nSpecie - 1;

        Info << nl << "yEquation.H: mixtureCheck=true. Number of species to solve: " << last << nl << endl;
        
        // Keep running sum for closure
        volScalarField sumY
        (
            IOobject
            (
                "sumY",
                runTime.timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh,
            dimensionedScalar(dimless, 0.0)
        );

        for (label i = 0; i < last; ++i)
        {
            volScalarField& Yi = Y[i];

            // Conservative update: ddt(rho*Yi) + div(phiRhoYi[i]) = 0
            // (phiRhoYi is already the mass flux of species i)
            solve
            (
                fvm::ddt(rho, Yi) + fvc::div(phiRhoYi[i])
            );

            // Basic bounding (important near shocks/contacts)
            Yi.max(0.0);
            Yi.min(1.0);

            sumY += Yi;
        }

        // Close last species and bound
        volScalarField& Ylast = Y[last];
        Ylast = scalar(1) - sumY;
        Ylast.max(0.0);
        Ylast.min(1.0);

        // Optional: renormalize so sum(Y)=1 exactly (after bounding)
        // This helps prevent drift when bounding clips values.
        volScalarField denom(sumY + Ylast);
        denom.max(SMALL);

        for (label i = 0; i < nSpecie; ++i)
        {
            Y[i] /= denom;
            Y[i].correctBoundaryConditions();
        }

        thermo2T.correct();
    }
}
