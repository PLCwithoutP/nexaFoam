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
    This file is for determining timestep and whether adjustable timestep is available.
    This file will be included to nexaFoam.C

\*---------------------------------------------------------------------------*/

// This block calculates the mean and maximum wave speed based Courant Numbers.

// #include "centralCourantNo.H"
if (mesh.nInternalFaces())
{
    scalarField sumAmaxSf(fvc::surfaceSum(max_a)().primitiveField());

    CoNum = 0.5*gMax(sumAmaxSf/mesh.V().field())*runTime.deltaTValue();

    meanCoNum =
        0.5*(gSum(sumAmaxSf)/gSum(mesh.V().field()))*runTime.deltaTValue();
}

Info<< "Mean and max Courant Numbers = "
    << meanCoNum << " " << CoNum << endl;


if (LTS) // LTS is created with time, I think (createTime.H)
{
    //#include "setRDeltaT.H"
    volScalarField& rDeltaT = trDeltaT.ref();

    scalar rDeltaTSmoothingCoeff
    (
        runTime.controlDict().getOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            0.02
        )
    );

    // Set the reciprocal time-step from the local Courant number
    rDeltaT.ref() = max
    (
        1/dimensionedScalar("maxDeltaT", dimTime, maxDeltaT),
        fvc::surfaceSum(max_a)()()
       /((2*maxCo)*mesh.V())
    );

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);

    Info << "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    ++runTime;
}

// #include "readTimeControls.H"
adjustTimeStep =
    runTime.controlDict().getOrDefault("adjustTimeStep", false);

maxCo =
    runTime.controlDict().getOrDefault<scalar>("maxCo", 1);

maxDeltaT =
    runTime.controlDict().getOrDefault<scalar>("maxDeltaT", GREAT);

if (!LTS)
{
    //#include "setDeltaT.H"
    if (adjustTimeStep)
    {
        scalar maxDeltaTFact = maxCo/(CoNum + SMALL);
        scalar deltaTFact = min(min(maxDeltaTFact, 1.0 + 0.1*maxDeltaTFact), 1.2);

        runTime.setDeltaT
        (
            min
            (
                deltaTFact*runTime.deltaTValue(),
                maxDeltaT
            )
        );

        Info << "deltaT = " <<  runTime.deltaTValue() << endl;
    }
    ++runTime;

    // Do any mesh changes
    mesh.update();
}