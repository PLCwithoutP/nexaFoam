/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "heNe2TReactionThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculate
(
    const volScalarField& p,
    volScalarField& TTR,
    volScalarField& TVib,
    volScalarField& h,
    volScalarField& eT,
    volScalarField& eR,
    volScalarField& eVib,
    volScalarField& psi,
    volScalarField& mu,
    volScalarField& alphaTR,
    const bool doOldTimes
)
{
    // Note: update oldTimes before current time so that if TTR.oldTime() is
    // created from TTR, it starts from the unconverted TTR
    if (doOldTimes && (p.nOldTimes() || TTR.nOldTimes() || TVib.nOldTimes()))
    {
        calculate
        (
            p.oldTime(),
            TTR.oldTime(),
            TVib.oldTime(),
            h.oldTime(),
            eT.oldTime(),
            eR.oldTime(),
            eVib.oldTime(),
            psi.oldTime(),
            mu.oldTime(),
            alphaTR.oldTime(),
            true
        );
    }

    // Mixing Rule definition for mixture
    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();
    WilkeMR<mixingMixtureType> wilkeMix(mix);

    PtrList<volScalarField>& Y_species = mix.Y();
    const speciesTable& spec  = mix.species();

    const scalar thetaVib_ = this->cellMixture(0).ThetaVib(); 
    const scalarField& hCells = h.primitiveField();
    const scalarField& eVibCells = eVib.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TTRCells = TTR.primitiveFieldRef();
    scalarField& TVibCells = TVib.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaTRCells = alphaTR.primitiveFieldRef();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);                  
        
        if (this->updateTTR())
        {
            TTRCells[celli] = cellMixture_.TH_TR
            (
                hCells[celli],
                pCells[celli],
                TTRCells[celli]
            );
        }

        if (this->updateTVib())
        {
            TVibCells[celli] = cellMixture_.TE_Vib
            (
                eVibCells[celli],
                pCells[celli],
                TTRCells[celli],
                TVibCells[celli]
            );
        }
        psiCells[celli] = cellMixture_.psi(pCells[celli], TTRCells[celli]);
        muCells[celli] = wilkeMix.muCell(celli, pCells[celli], TTRCells[celli]);
        alphaTRCells[celli] = wilkeMix.alphaTRCell(celli, pCells[celli], TTRCells[celli]);
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    volScalarField::Boundary& TTRBf = TTR.boundaryFieldRef();
    volScalarField::Boundary& TVibBf = TVib.boundaryFieldRef();
    volScalarField::Boundary& psiBf = psi.boundaryFieldRef();
    volScalarField::Boundary& hBf = h.boundaryFieldRef();
    volScalarField::Boundary& eTBf = eT.boundaryFieldRef();
    volScalarField::Boundary& eRBf = eR.boundaryFieldRef();
    volScalarField::Boundary& eVibBf = eVib.boundaryFieldRef();
    volScalarField::Boundary& muBf = mu.boundaryFieldRef();
    volScalarField::Boundary& alphaBf = alphaTR.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        fvPatchScalarField& pTTR = TTRBf[patchi];
        fvPatchScalarField& pTVib = TVibBf[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];
        fvPatchScalarField& ph = hBf[patchi];
        fvPatchScalarField& pesT = eTBf[patchi];
        fvPatchScalarField& pesR = eRBf[patchi];
        fvPatchScalarField& pesVib = eVibBf[patchi];
        fvPatchScalarField& pmu = muBf[patchi];
        fvPatchScalarField& palphaTR = alphaBf[patchi];

        const fvPatch& patch = pp.patch();
        const labelList& fc = patch.faceCells();
        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const label celli = fc[facei];
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = cellMixture_.H(pp[facei], pTTR[facei]);
                pesT[facei] = cellMixture_.ET(pp[facei], pTTR[facei]);
                pesR[facei] = cellMixture_.ER(pp[facei], pTTR[facei]);
                pesVib[facei] = cellMixture_.EV(pp[facei], pTTR[facei], pTVib[facei], thetaVib_);

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);
                pmu[facei] = wilkeMix.muCell(celli, pp[facei], pTTR[facei]);
                palphaTR[facei] = wilkeMix.alphaTRCell(celli, pp[facei], pTTR[facei]);
            }
        }
        else
        {
            forAll(pTTR, facei)
            {
                const label celli = fc[facei];
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateTTR())
                {
                    pTTR[facei] = cellMixture_.TH_TR(ph[facei], pp[facei], pTTR[facei]);
                }

                if (this->updateTVib())
                {
                    pTVib[facei] = cellMixture_.TE_Vib(pesVib[facei], pp[facei], pTTR[facei], pTVib[facei]);
                }

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);
                pmu[facei] = wilkeMix.muCell(celli, pp[facei], pTTR[facei]);
                palphaTR[facei] = wilkeMix.alphaTRCell(celli, pp[facei], pTTR[facei]);
            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateTEnergy
(
    const volScalarField& p,
    const volScalarField& TTR,
    volScalarField& eT
)
{

    scalarField& eTCells = eT.primitiveFieldRef();

    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);
        
        eTCells[celli] = cellMixture_.ET
        (
            pCells[celli],
            TTRCells[celli]
        );

    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TTRBf = TTR.boundaryField();
    volScalarField::Boundary& eTBf = eT.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pTTR = TTRBf[patchi];
        fvPatchScalarField& pesT = eTBf[patchi];

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                pesT[facei] = cellMixture_.ET(pp[facei], pTTR[facei]);

            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateREnergy
(
    const volScalarField& p,
    const volScalarField& TTR,
    volScalarField& eR
)
{

    scalarField& eRCells = eR.primitiveFieldRef();

    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        eRCells[celli] = cellMixture_.ER
        (
            pCells[celli],
            TTRCells[celli]
        );
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TTRBf = TTR.boundaryField();
    volScalarField::Boundary& eRBf = eR.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pTTR = TTRBf[patchi];
        fvPatchScalarField& pesR = eRBf[patchi];

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                pesR[facei] = cellMixture_.ER(pp[facei], pTTR[facei]);

            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateVibEnergy
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    volScalarField& eVib
)
{
    const scalar thetaVib_ = this->cellMixture(0).ThetaVib(); 

    scalarField& eVibCells = eVib.primitiveFieldRef();

    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        eVibCells[celli] = cellMixture_.EV
        (
            pCells[celli],
            TTRCells[celli],
            TVibCells[celli],
            thetaVib_
        );
    }

    const volScalarField::Boundary& pBf = p.boundaryField();
    const volScalarField::Boundary& TTRBf = TTR.boundaryField();
    const volScalarField::Boundary& TVibBf = TVib.boundaryField();
    volScalarField::Boundary& eVibBf = eVib.boundaryFieldRef();

    forAll(pBf, patchi)
    {
        const fvPatchScalarField& pp = pBf[patchi];
        const fvPatchScalarField& pTTR = TTRBf[patchi];
        const fvPatchScalarField& pTVib = TVibBf[patchi];
        fvPatchScalarField& pesVib = eVibBf[patchi];

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                pesVib[facei] = cellMixture_.EV(pp[facei], pTTR[facei], pTVib[facei], thetaVib_);

            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::wilkeKappaAverage
() const
{
    // Mixing Rule definition for mixture
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();
    WilkeMR<mixingMixtureType> wilkeMix(mix);

    const scalarField& TTRCells = this->TTR_.primitiveField();
    const scalarField& pCells = this->p_.primitiveField();

    const fvMesh& mesh = this->TTR_.mesh();
    const scalarField& V = mesh.V();
    const scalar Vtot = gSum(V);

    scalar sumKappaV = 0.0;

    forAll(TTRCells, celli)
    {
        const scalar kappaCell =
            wilkeMix.kappaTRCell(celli, pCells[celli], TTRCells[celli]);

        sumKappaV += kappaCell * V[celli];
    }

    return sumKappaV / Vtot;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::wilkeMuAverage
() const
{
    // Mixing Rule definition for mixture
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();
    WilkeMR<mixingMixtureType> wilkeMix(mix);

    const scalarField& TTRCells = this->TTR_.primitiveField();
    const scalarField& pCells = this->p_.primitiveField();

    const fvMesh& mesh = this->TTR_.mesh();
    const scalarField& V = mesh.V();
    const scalar Vtot = gSum(V);

    scalar sumMuV = 0.0;

    forAll(TTRCells, celli)
    {
        
        const scalar muCell =
            wilkeMix.muCell(celli, pCells[celli], TTRCells[celli]);

        sumMuV += muCell * V[celli];
    }

    return sumMuV / Vtot;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateRhoMixture
(
    const scalar pCell,
    const scalar TTRCell,
    const label  celli
) const
{
    const typename MixtureType::thermoType& mixCell =
        this->cellMixture(celli);

    return mixCell.rho(pCell, TTRCell);
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::calculateCpTRMixture
(
    const scalar pCell,
    const scalar TTRCell,
    const label  celli
) const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& comp = this->composition();

    scalar cpMix = 0.0;

    // PtrList of Y-fields, one per species
    const auto& Y = comp.Y();
    const label nSpec = Y.size();

    for (label i = 0; i < nSpec; ++i)
    {
        const scalar Yi   = Y[i][celli];   // mass fraction of species i in this cell

        const scalar cp_i = comp.CpTR(i, pCell, TTRCell);

        cpMix += Yi*cp_i;
    }

    return cpMix;   // mixture cp_TR [J/(kg·K)] if cp_i is mass-based
}


template<class BasicNe2TThermo, class MixtureType>
Foam::volScalarField
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::fickDiffusionCoeff() const
{
    using mixingMixtureType = const typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mixingComp = this->composition();
    WilkeMR<mixingMixtureType> wilkeMix(mixingComp);

    const fvMesh& mesh = this->TTR_.mesh();

    volScalarField Deff
    (
        IOobject
        (
            "Deff",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimArea/dimTime, 0.0)   // [m^2/s]
    );

    FickDM<mixingMixtureType> fickDM(mixingComp);

    forAll(Deff, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);
        const scalar pCell      = this->p_[celli];
        const scalar TTRCell    = this->TTR_[celli];
        const scalar CpTRCell =
            this->calculateCpTRMixture(pCell, TTRCell, celli);
        const scalar kappaTRCell =
            wilkeMix.kappaTRCell(celli, pCell, TTRCell);
        const scalar rhoCell = 
            this->calculateRhoMixture(pCell, TTRCell, celli);
        Info << "CpTR cell is : " << CpTRCell << " in cell " << celli << nl;
        Info << "kappaTR cell is : " << kappaTRCell << " in cell " << celli << nl;
        Info << "rho cell is : " << rhoCell << " in cell " << celli << nl;
        Deff[celli] = fickDM.DeffCell
        (
            celli,
            pCell,
            TTRCell,
            rhoCell,
            kappaTRCell,
            CpTRCell
        );
        Info << "Deff is : " << Deff[celli] << " in cell " << celli << nl;
    }

    return Deff;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::heNe2TReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    he2TThermo<BasicNe2TThermo, MixtureType>(mesh, phaseName)
{
    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->h_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}


template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::heNe2TReactionThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    he2TThermo<BasicNe2TThermo, MixtureType>(mesh, phaseName, dictionaryName)
{
    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->h_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        true                    // Create old time fields
    );
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::~heNe2TReactionThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correct()
{
    DebugInFunction << endl;

    calculate
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->h_,
        this->eT_,
        this->eR_,
        this->eVib_,
        this->psi_,
        this->mu_,
        this->alpha_,
        false           // No need to update old times
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctTEnergy()
{
    DebugInFunction << endl;

    calculateTEnergy
    (
        this->p_,
        this->TTR_,
        this->eT_
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctREnergy()
{
    DebugInFunction << endl;

    calculateREnergy
    (
        this->p_,
        this->TTR_,
        this->eR_
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibEnergy()
{
    DebugInFunction << endl;

    calculateVibEnergy
    (
        this->p_,
        this->TTR_,
        this->TVib_,
        this->eVib_
    );

    DebugInFunction << "Finished" << endl;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibSource()
{
    DebugInFunction << endl;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    static VTEnergySource<mixingMixtureType, mixingRuleType> Q_vt_source(mix, wilkeMix);

    PtrList<volScalarField>& Qlist =
        Q_vt_source.correctVibSource(this->p_, this->TTR_, this->TVib_);

    DebugInFunction << "Finished" << endl;

    return Qlist;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVibVibSource()
{
    DebugInFunction << endl;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    static VVEnergySource<mixingMixtureType, mixingRuleType> Q_vv_source(mix, wilkeMix);

    PtrList<volScalarField>& QVVlist =
        Q_vv_source.correctVibVibSource(this->p_, this->TTR_, this->TVib_);

    DebugInFunction << "Finished" << endl;

    return QVVlist;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TReactionThermo<BasicNe2TThermo, MixtureType>::correctVTRelaxationTime()
{
    DebugInFunction << endl;

    using mixingMixtureType = typename MixtureType::basicSpecie2TMixture;
    mixingMixtureType& mix = this->composition();

    using mixingRuleType = Foam::WilkeMR<mixingMixtureType>;

    static mixingRuleType wilkeMix(mix);
    static VTEnergySource<mixingMixtureType, mixingRuleType> Q_vt_source(mix, wilkeMix);

    PtrList<volScalarField>& TauList =
        Q_vt_source.correctVTRelaxationTime(this->p_, this->TTR_);

    DebugInFunction << "Finished" << endl;

    return TauList;
}

// ************************************************************************* //
