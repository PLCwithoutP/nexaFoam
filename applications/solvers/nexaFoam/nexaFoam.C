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
#include "fixed2TRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "FluxCalculators/kurganov.H"
#include "DiffusionModels/FickModel.H"
#include "BasicNeChemistryModel.H"

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
    #include "createTimeControls.H"


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    const dimensionedScalar v_zero(dimVolume/dimTime, Zero);

    // Flux Calculator
    autoPtr<Foam::kurganovFluxCalculator> kurganovFluxPtr(nullptr);

    if (!use2T)
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                hTR,
                TTR,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }
    else if (use2T && !mixtureCheck)
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                hTR,
                TTR,
                eV,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }
    else
    {
        kurganovFluxPtr.reset
        (
            new Foam::kurganovFluxCalculator
            (
                mesh,
                thermo2T,
                rho,
                rhoU,
                U,
                psi,
                hTR,
                TTR,
                eVibSpecies,
                pos,
                neg,
                phi,
                fluxScheme,
                v_zero,
                mixtureCheck,
                use2T
            )
        );
    }

    Foam::kurganovFluxCalculator& kurganovFlux = kurganovFluxPtr();

    surfaceVectorField& phiUp = kurganovFlux.phiUp();
    surfaceScalarField& phiEp = kurganovFlux.phiEp();
    surfaceScalarField& sigmaDotU = kurganovFlux.sigmaDotU();
    surfaceScalarField& max_a = kurganovFlux.maxA();
    surfaceScalarField& phiEv = kurganovFlux.phiEv();
    const PtrList<surfaceScalarField>& phiRhoEvibYi = kurganovFlux.phiRhoEvibYi();

    // Courant numbers used to adjust the time-step
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        //#include "runTimeLogging.H"

        volScalarField muLam("muLam", thermo2T.mu());
        volTensorField tauMC("tauMC", muLam*dev2(Foam::T(fvc::grad(U))));

        // Effective diffusion coefficient must be updated in every timestep
        Deff = thermo2T.fickDiffusionCoeff();
        Deff.correctBoundaryConditions();

        kurganovFlux.invoke(muLam, tauMC);

        #include "decideDeltaT.H"

        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "Equations/continuityEquation.H"

        #include "Equations/yEquation.H"

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
              - fvm::laplacian(muLam, U)
              - fvc::div(tauMC)
            );
            rhoU = rho*U;
        }

        #include "Equations/totalEnergyEquation.H"
        
        #include "Equations/vibEnergyEquation.H"

        #include "Updates/updateTTRTemperature.H"

        // Sanity check
        if (mixtureCheck)
        {
            volScalarField drho(rhoCont - rho);
            const scalar maxErr = gMax(mag(drho)().primitiveField());
            Info<< "Sanity check: max(|rhoCont - rhoFromSpecies|) = "
                << maxErr << nl << endl;
        }

        #include "Updates/updateVibTemperature.H"

        #include "Equations/applyChemistry.H"

        runTime.write();

        runTime.printExecutionTime(Info);

    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
