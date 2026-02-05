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
    Mass diffusion flux calculation of species.
    J_s,j is calculated in this snippet

\*---------------------------------------------------------------------------*/


// Mass Diffusion Flux at Face-centers : for equations
tmp<surfaceScalarField> fickFluxFace
(
    const volScalarField& rho,
    const volScalarField& Deff,
    const volScalarField& Ys
)
{
    const fvMesh& mesh = Ys.mesh();

    // Interpolate coefficient to faces
    tmp<surfaceScalarField> tRhoDeff_f = fvc::interpolate(rho*Deff);
    const surfaceScalarField& rhoDeff_f = tRhoDeff_f();

    // Normal gradient of Ys on faces (1/m)
    tmp<surfaceScalarField> tSnGradYs = fvc::snGrad(Ys);
    const surfaceScalarField& snGradYs = tSnGradYs();

    // J·n on faces (kg / (m^2 s))
    tmp<surfaceScalarField> tJs_n
    (
        new surfaceScalarField
        (
            IOobject
            (
                "Js_n_" + Ys.name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            - rhoDeff_f * snGradYs
        )
    );

    return tJs_n;
}

tmp<surfaceScalarField> tJs_n = fickFluxFace(rho, Deff_s, Ys);
surfaceScalarField phiJ = tJs_n()*mesh.magSf();
