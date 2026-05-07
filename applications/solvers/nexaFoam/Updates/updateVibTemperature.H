//Info<< "Entering updateTVibTemperature" << nl << endl;
// Pure mixture 2T case
if (use2T && !mixtureCheck)
{
    thermo2T.correctTVib();
    TVib.correctBoundaryConditions();

    eV.correctBoundaryConditions();
    rhoEv = rho*eV;
    rhoEv.correctBoundaryConditions();
}

// Multi-specie 2T case
if (use2T && mixtureCheck)
{
    const scalarField& pCells   = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    scalarField& TVibMixCells   = TVib.primitiveFieldRef();
    scalarField& eVCells        = eV.primitiveFieldRef();

    forAll(Y_species, speciei)
    {
        if (!vibSpecieActive[speciei])
        {
            TVibSpecies[speciei] =
                dimensionedScalar("zero", dimTemperature, 0.0);
            TVibSpecies[speciei].correctBoundaryConditions();
            continue;
        }

        scalarField& TVs = TVibSpecies[speciei].primitiveFieldRef();
        const scalarField& eVs = eVibSpecies[speciei].primitiveField();

        forAll(TVs, celli)
        {
            const scalar eVsClamped = max(eVs[celli], scalar(0));

            if (eVsClamped < SMALL)
            {
                // SHO ground state: no vibrational excitation, TVib → 0
                TVs[celli] = scalar(0);
                continue;
            }

            TVs[celli] =
                mix.TEVib
                (
                    speciei,
                    eVsClamped,
                    pCells[celli],
                    TTRCells[celli],
                    TVs[celli]
                );
        }

        const volScalarField::Boundary& pBf   = p.boundaryField();
        const volScalarField::Boundary& TTRBf = TTR.boundaryField();
        const volScalarField::Boundary& eVBf  = eVibSpecies[speciei].boundaryField();
        volScalarField::Boundary& TVBf        = TVibSpecies[speciei].boundaryFieldRef();

        forAll(pBf, patchi)
        {
            const fvPatchScalarField& pp   = pBf[patchi];
            const fvPatchScalarField& pTTR = TTRBf[patchi];
            const fvPatchScalarField& peV  = eVBf[patchi];
            fvPatchScalarField& pTV        = TVBf[patchi];

            forAll(pp, facei)
            {
                const scalar eVsClamped = max(peV[facei], scalar(0));

                if (eVsClamped < SMALL)
                {
                    pTV[facei] = scalar(0);
                    continue;
                }

                pTV[facei] =
                    mix.TEVib
                    (
                        speciei,
                        eVsClamped,
                        pp[facei],
                        pTTR[facei],
                        pTV[facei]
                    );
            }
        }

        TVibSpecies[speciei].correctBoundaryConditions();
    }

    // Reconstruct mixture TVib only as a compatibility field for the current
    // thermo/source infrastructure, which still expects one TVib field.
    forAll(TVibMixCells, celli)
    {
        scalar sumY = 0.0;
        scalar sumTV = 0.0;

        forAll(Y_species, speciei)
        {
            if (!vibSpecieActive[speciei])
            {
                continue;
            }

            const scalar Yi = max(Y_species[speciei][celli], scalar(0));
            sumY += Yi;
            sumTV += Yi*TVibSpecies[speciei][celli];
        }

        TVibMixCells[celli] = (sumY > SMALL ? sumTV/sumY : 0.0);
    }

    forAll(TVib.boundaryField(), patchi)
    {
        fvPatchScalarField& pTVmix = TVib.boundaryFieldRef()[patchi];

        forAll(pTVmix, facei)
        {
            scalar sumY = 0.0;
            scalar sumTV = 0.0;

            forAll(Y_species, speciei)
            {
                if (!vibSpecieActive[speciei])
                {
                    continue;
                }

                const scalar Yi =
                    max(Y_species[speciei].boundaryField()[patchi][facei], scalar(0));

                sumY += Yi;
                sumTV += Yi*TVibSpecies[speciei].boundaryField()[patchi][facei];
            }

            pTVmix[facei] = (sumY > SMALL ? sumTV/sumY : 0.0);
        }
    }

    TVib.correctBoundaryConditions();

    // ── Mutation++ path: add atomic electronic energy diagnostically ────────────
    // Atoms have no VT source in the FVM equation (vibSpecieActive=false).
    // Their electronic energy is passively determined by T_ve, not independently evolved.
    // Update eVibSpecies[atom] = e_el(atom, TVib) and correct rhoEv to include it.
    if (use2T && mixtureCheck && thermo2T.hasElectronicVibCoupling())
    {
        forAll(Y_species, speciei)
        {
            if (mix.isSpecieMolecular(speciei)) continue;

            scalarField& eVsAtom  = eVibSpecies[speciei].primitiveFieldRef();
            scalarField& rYeVAtom = rhoYEvibSpecies[speciei].primitiveFieldRef();
            const scalarField& TVibCells = TVib.primitiveField();

            forAll(TVibCells, celli)
            {
                eVsAtom[celli] = mix.EsVib
                (
                    speciei, p[celli], TTR[celli], TVibCells[celli], 0.0
                );
                rYeVAtom[celli] = rhoYi[speciei][celli] * eVsAtom[celli];
            }
            eVibSpecies[speciei].correctBoundaryConditions();
            rhoYEvibSpecies[speciei].correctBoundaryConditions();
        }

        // Recompute total eV including atomic electronic contribution
        eV = dimensionedScalar("zero", eV.dimensions(), 0.0);
        forAll(Y_species, speciei)
        {
            if (vibSpecieActive[speciei] || !mix.isSpecieMolecular(speciei))
            {
                eV += rhoYEvibSpecies[speciei] / rho;
            }
        }
        eV.correctBoundaryConditions();
        rhoEv = rho * eV;
        rhoEv.correctBoundaryConditions();
    }
}
