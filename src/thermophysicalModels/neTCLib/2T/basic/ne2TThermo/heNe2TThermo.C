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

#include "heNe2TThermo.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::calculate
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
    volScalarField& alpha,
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
            alpha.oldTime(),
            true
        );
    }

    const bool use2T = this->twoTemperature();
    
    const scalar thetaVib_ = this->cellMixture(0).ThetaVib(); 
    const scalarField& hCells = h.primitiveField();
    const scalarField& eVibCells = eVib.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& TTRCells = TTR.primitiveFieldRef();
    scalarField& TVibCells = TVib.primitiveFieldRef();
    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& muCells = mu.primitiveFieldRef();
    scalarField& alphaCells = alpha.primitiveFieldRef();

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

        if (!use2T)
        {
            TVibCells[celli] = TTRCells[celli];
        }
        else if (this->updateTVib())
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
        muCells[celli] = cellMixture_.mu(TTRCells[celli]);
        alphaCells[celli] = cellMixture_.alphah(pCells[celli], TTRCells[celli]);
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
    volScalarField::Boundary& alphaBf = alpha.boundaryFieldRef();

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
        fvPatchScalarField& palpha = alphaBf[patchi];

        if (pTTR.fixesValue())
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                ph[facei] = cellMixture_.H(pp[facei], pTTR[facei]);
                pesT[facei] = cellMixture_.ET(pp[facei], pTTR[facei]);
                pesR[facei] = cellMixture_.ER(pp[facei], pTTR[facei]);

                if (!use2T)
                {
                    pTVib[facei] = pTTR[facei];
                }

                const scalar TVibUse = use2T ? pTVib[facei] : pTTR[facei];

                pesVib[facei] = cellMixture_.EV
                (
                    pp[facei],
                    pTTR[facei],
                    TVibUse,
                    thetaVib_
                );

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);
                pmu[facei] = cellMixture_.mu(pTTR[facei]);
                palpha[facei] = cellMixture_.alphah(pp[facei], pTTR[facei]);
            }
        }
        else
        {
            forAll(pTTR, facei)
            {
                const typename MixtureType::thermoType& cellMixture_ =
                    this->patchFaceMixture(patchi, facei);

                if (this->updateTTR())
                {
                    pTTR[facei] = cellMixture_.TH_TR
                    (
                        ph[facei], 
                        pp[facei], 
                        pTTR[facei]
                    );
                }

                if (!use2T)
                {
                    pTVib[facei] = pTTR[facei];
                }
                else if (this->updateTVib())
                {
                    pTVib[facei] = cellMixture_.TE_Vib
                    (
                        pesVib[facei],
                        pp[facei],
                        pTTR[facei],
                        pTVib[facei]
                    );
                }

                ppsi[facei] = cellMixture_.psi(pp[facei], pTTR[facei]);
                pmu[facei] = cellMixture_.mu(pTTR[facei]);
                palpha[facei] = cellMixture_.alphah(pp[facei], pTTR[facei]);
            }
        }
    }
}

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::calculateTEnergy
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
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::calculateREnergy
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
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::calculateVibEnergy
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    volScalarField& eVib
)
{
    const bool use2T = this->twoTemperature();
    const scalar thetaVib_ = this->cellMixture(0).ThetaVib(); 

    scalarField& eVibCells = eVib.primitiveFieldRef();

    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(TTRCells, celli)
    {
        const typename MixtureType::thermoType& cellMixture_ =
            this->cellMixture(celli);

        const scalar TVibUse = use2T ? TVibCells[celli] : TTRCells[celli];

        eVibCells[celli] = cellMixture_.EV
        (
            pCells[celli],
            TTRCells[celli],
            TVibUse,
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

                const scalar TVibUse = use2T ? pTVib[facei] : pTTR[facei];

                pesVib[facei] = cellMixture_.EV
                (
                    pp[facei],
                    pTTR[facei],
                    TVibUse,
                    thetaVib_
                );

            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::heNe2TThermo
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
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::heNe2TThermo
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
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::~heNe2TThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicNe2TThermo, class MixtureType>
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correct()
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
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctTEnergy()
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
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctREnergy()
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
void Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctVibEnergy()
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
Foam::scalar
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::wilkeKappaAverage
() const
{
    return 0;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::scalar
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::wilkeMuAverage
() const
{
    return 0;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>& Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctVibSource
()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList;
    }

    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>& Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctVibVibSource
()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList;
    }

    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>&
Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctCVSource
(
    const PtrList<volScalarField::Internal>& RR
)
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList;
    }

    // Non-reacting thermo: no C-V chemical source
    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::PtrList<Foam::volScalarField>& Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::correctVTRelaxationTime
()
{
    static Foam::PtrList<Foam::volScalarField> emptyList;

    if (!this->twoTemperature())
    {
        return emptyList;
    }

    return emptyList;
}

template<class BasicNe2TThermo, class MixtureType>
Foam::volScalarField Foam::heNe2TThermo<BasicNe2TThermo, MixtureType>::fickDiffusionCoeff
() const
{
    const fvMesh& mesh = this->TTR_.mesh();  

    return volScalarField
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
        dimensionedScalar("zero", dimArea/dimTime, 0.0)  // [m^2/s]
    );
}
// ************************************************************************* //
