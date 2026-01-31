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

Application
    nexaFoam

Group
    grpCompressibleSolvers

Description
    Density-based compressible flow solver based on
    central-upwind schemes of Kurganov and Tadmor with
    support for mesh-motion and topology changes.

\*---------------------------------------------------------------------------*/


#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "ne2TReactionThermo.H"
#include "turbulentFluidTThermoModel.H"
#include "fixed2TRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "logDebugCreate.H"
    argList::addNote
    (
        "Density-based compressible flow solver based on"
        " central-upwind schemes of Kurganov and Tadmor with"
        " support for mesh-motion and topology changes."
    );

    #define NO_CONTROL
    #include "postProcess.H"
    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createTimeControls.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //#include "runTimeLogging.H"
        //Info << "Effective diffusion coefficient of mixture is : \n" << (thermo2T.fickDiffusionCoeff())().average().value() << endl;

        volScalarField muEff("muEff", turbulence->muEff());
        volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));

        #include "FluxCalculators/kurganov.H"

        #include "decideDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Solve density
        solve(fvm::ddt(rho) + fvc::div(phi));

        // --- Solve momentum
        solve(fvm::ddt(rhoU) + fvc::div(phiUp));

        U.ref() =
            rhoU()
           /rho();
        U.correctBoundaryConditions();
        rhoU.boundaryFieldRef() == rho.boundaryField()*U.boundaryField();

        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, U) - fvc::ddt(rho, U)
              - fvm::laplacian(muEff, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        solve
        (
            fvm::ddt(rhoE)
          + fvc::div(phiEp)
          - fvc::div(sigmaDotU)
        );

        hTR = rhoE/rho - 0.5*magSqr(U);
        hTR.correctBoundaryConditions();
        eT.correctBoundaryConditions();
        if (vibrationalCheck)
        {
            eR.correctBoundaryConditions();
            eV.correctBoundaryConditions();
        }
        thermo2T.correct();
        thermo2T.correctTEnergy();
        if (vibrationalCheck)
        {
            thermo2T.correctREnergy();
            thermo2T.correctVibEnergy();
        }
        rhoE.boundaryFieldRef() ==
            rho.boundaryField()*
            (
                hTR.boundaryField() + 0.5*magSqr(U.boundaryField())
            );
        if (!inviscid)
        {
            solve
            (
                fvm::ddt(rho, hTR) - fvc::ddt(rho, hTR)
              - fvm::laplacian(turbulence->alphaEff(), hTR)
            );
            thermo2T.correct();
            thermo2T.correctTEnergy();
            if (vibrationalCheck)
            {
                thermo2T.correctREnergy();
                thermo2T.correctVibEnergy();
            }
            rhoE = rho*(hTR + 0.5*magSqr(U));        
        }
        p.ref() =
            rho()
           /psi();
        p.correctBoundaryConditions();
        rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

        turbulence->correct();

        runTime.write();

        runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
