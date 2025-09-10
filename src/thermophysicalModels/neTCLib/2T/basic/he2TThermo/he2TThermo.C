/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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

#include "he2TThermo.H"
#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
void Foam::he2TThermo<Basic2TThermo, MixtureType>::
heBoundaryCorrection(volScalarField& h)
{
    volScalarField::Boundary& hBf = h.boundaryFieldRef();

    forAll(hBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(hBf[patchi]).gradient()
                = hBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(hBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(hBf[patchi]).refGrad()
                = hBf[patchi].fvPatchField::snGrad();
        }
    }
}


template<class Basic2TThermo, class MixtureType>
void Foam::he2TThermo<Basic2TThermo, MixtureType>::init
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib,
    volScalarField& he
)
{
    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p.primitiveField();
    const scalarField& TTRCells = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        heBf[patchi] == this->he
        (
            p.boundaryField()[patchi],
            TTR.boundaryField()[patchi],
            patchi
        );

        heBf[patchi].useImplicit(TTR.boundaryField()[patchi].useImplicit());
    }

    this->heBoundaryCorrection(he);

    // Note: TTR does not have oldTime
    if (p.nOldTimes())
    {
        init(p.oldTime(), TTR.oldTime(), TVib.oldTime(), he.oldTime());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const word& phaseName
)
:
    Basic2TThermo(mesh, phaseName),
    MixtureType(*this, mesh, phaseName),

    he_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init(this->p_, this->TTR_, this->TVib_, he_);
}


template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const dictionary& dict,
    const word& phaseName
)
:
    Basic2TThermo(mesh, dict, phaseName),
    MixtureType(*this, mesh, phaseName),

    he_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init(this->p_, this->TTR_, this->TVib_, he_);
}


template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::he2TThermo
(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName
)
:
    Basic2TThermo(mesh, phaseName, dictionaryName),
    MixtureType(*this, mesh, phaseName),

    he_
    (
        IOobject
        (
            Basic2TThermo::phasePropertyName
            (
                MixtureType::thermoType::heName()
            ),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimEnergy/dimMass,
        this->heBoundaryTypes(),
        this->heBoundaryBaseTypes()
    )
{
    init(this->p_, this->TTR_, this->TVib_, he_);
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::he2TThermo<Basic2TThermo, MixtureType>::~he2TThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::he
(
    const volScalarField& p,
    const volScalarField& TTR
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto the = volScalarField::New
    (
        "he",
        IOobject::NO_REGISTER,
        mesh,
        he_.dimensions()
    );
    auto& he = the.ref();

    scalarField& heCells = he.primitiveFieldRef();
    const scalarField& pCells = p;
    const scalarField& TTRCells = TTR;

    forAll(heCells, celli)
    {
        heCells[celli] =
            this->cellMixture(celli).HE(pCells[celli], TTRCells[celli]);
    }

    volScalarField::Boundary& heBf = he.boundaryFieldRef();

    forAll(heBf, patchi)
    {
        scalarField& hep = heBf[patchi];
        const scalarField& pp = p.boundaryField()[patchi];
        const scalarField& TTRp = TTR.boundaryField()[patchi];

        forAll(hep, facei)
        {
            hep[facei] =
                this->patchFaceMixture(patchi, facei).HE(pp[facei], TTRp[facei]);
        }
    }

    return the;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto the = tmp<scalarField>::New(TTR.size());
    auto& he = the.ref();

    forAll(TTR, celli)
    {
        he[celli] = this->cellMixture(cells[celli]).HE(p[celli], TTR[celli]);
    }

    return the;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::he
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto the = tmp<scalarField>::New(TTR.size());
    auto& he = the.ref();

    forAll(TTR, facei)
    {
        he[facei] =
            this->patchFaceMixture(patchi, facei).HE(p[facei], TTR[facei]);
    }

    return the;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::hc() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto thc = volScalarField::New
    (
        "hc",
        IOobject::NO_REGISTER,
        mesh,
        he_.dimensions()
    );
    auto& hcf = thc.ref();

    scalarField& hcCells = hcf.primitiveFieldRef();

    forAll(hcCells, celli)
    {
        hcCells[celli] = this->cellMixture(celli).Hc();
    }

    volScalarField::Boundary& hcfBf = hcf.boundaryFieldRef();

    forAll(hcfBf, patchi)
    {
        scalarField& hcp = hcfBf[patchi];

        forAll(hcp, facei)
        {
            hcp[facei] = this->patchFaceMixture(patchi, facei).Hc();
        }
    }

    return thc;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpTR = tmp<scalarField>::New(TTR.size());
    auto& cpTR = tCpTR.ref();

    forAll(TTR, facei)
    {
        cpTR[facei] =
            this->patchFaceMixture(patchi, facei).CpTR(p[facei], TTR[facei]);
    }

    return tCpTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tCpTR = tmp<scalarField>::New(TTR.size());
    auto& CpTR = tCpTR.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        CpTR[i] = this->cellMixture(celli).CpTR(p[i], TTR[i]);
    }

    return tCpTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpTR = volScalarField::New
    (
        "CpTR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cpTR = tCpTR.ref();

    forAll(this->TTR_, celli)
    {
        cpTR[celli] =
            this->cellMixture(celli).CpTR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cpTRBf = cpTR.boundaryFieldRef();

    forAll(cpTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpTR = cpTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpTR[facei] =
                this->patchFaceMixture(patchi, facei).CpTR(pp[facei], pTTR[facei]);
        }
    }

    return tCpTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvT
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCvT = tmp<scalarField>::New(TTR.size());
    auto& cvT = tCvT.ref();

    forAll(TTR, facei)
    {
        cvT[facei] =
            this->patchFaceMixture(patchi, facei).CvT(p[facei], TTR[facei]);
    }

    return tCvT;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCvR = tmp<scalarField>::New(TTR.size());
    auto& cvR = tCvR.ref();

    forAll(TTR, facei)
    {
        cvR[facei] =
            this->patchFaceMixture(patchi, facei).CvR(p[facei], TTR[facei]);
    }

    return tCvR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvVib
(
    const scalarField& p,
    const scalarField& TTR,
    const scalarField& TVib,
    const scalar ThetaVib,
    const label patchi
) const
{
    auto tCvVib = tmp<scalarField>::New(TVib.size());
    auto& cvVib = tCvVib.ref();

    forAll(TVib, facei)
    {
        cvVib[facei] =
            this->patchFaceMixture(patchi, facei).CvVib(p[facei], TTR[facei], TVib[facei], ThetaVib);
    }

    return tCvVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::rhoEoS
(
    const scalarField& p,
    const scalarField& TTR,
    const labelList& cells
) const
{
    auto tRho = tmp<scalarField>::New(TTR.size());
    auto& rho = tRho.ref();

    forAll(cells, i)
    {
        const label celli = cells[i];
        rho[i] = this->cellMixture(celli).rho(p[i], TTR[i]);
    }

    return tRho;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvT() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCvT = volScalarField::New
    (
        "CvT",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvT = tCvT.ref();

    forAll(this->TTR_, celli)
    {
        cvT[celli] =
            this->cellMixture(celli).CvT(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cvTBf = cvT.boundaryFieldRef();

    forAll(cvTBf, patchi)
    {
        cvTBf[patchi] = CvT
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        );
    }

    return tCvT;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCvR = volScalarField::New
    (
        "CvR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvR = tCvR.ref();

    forAll(this->TTR_, celli)
    {
        cvR[celli] =
            this->cellMixture(celli).CvR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& cvRBf = cvR.boundaryFieldRef();

    forAll(cvRBf, patchi)
    {
        cvRBf[patchi] = CvR
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        );
    }

    return tCvR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CvVib() const
{
    const fvMesh& mesh = this->TVib_.mesh();

    auto tCvVib = volScalarField::New
    (
        "CvVib",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& cvVib = tCvVib.ref();

    forAll(this->TVib_, celli)
    {
        cvVib[celli] =
            this->cellMixture(celli).CvVib(this->p_[celli], this->TTR_[celli] ,this->TVib_[celli], this->ThetaVib_);
    }

    volScalarField::Boundary& cvVibBf = cvVib.boundaryFieldRef();

    forAll(cvVibBf, patchi)
    {
        cvVibBf[patchi] = CvVib
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            this->TVib_.boundaryField()[patchi],
            this->ThetaVib_,
            patchi
        );
    }

    return tCvVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::gamma
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tgamma = tmp<scalarField>::New(TTR.size());
    auto& gamma = tgamma.ref();

    forAll(TTR, facei)
    {
        gamma[facei] =
            this->patchFaceMixture(patchi, facei).gamma(p[facei], TTR[facei]);
    }

    return tgamma;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::gamma() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tgamma = volScalarField::New
    (
        "gamma",
        IOobject::NO_REGISTER,
        mesh,
        dimless
    );
    auto& gamma = tgamma.ref();

    forAll(this->TTR_, celli)
    {
        gamma[celli] =
            this->cellMixture(celli).gamma(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& gammaBf = gamma.boundaryFieldRef();

    forAll(gammaBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pgamma = gammaBf[patchi];

        forAll(pTTR, facei)
        {
            pgamma[facei] = this->patchFaceMixture(patchi, facei).gamma
            (
                pp[facei],
                pTTR[facei]
            );
        }
    }

    return tgamma;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpvTR = tmp<scalarField>::New(TTR.size());
    auto& CpvTR = tCpvTR.ref();

    forAll(TTR, facei)
    {
        CpvTR[facei] =
            this->patchFaceMixture(patchi, facei).CpvTR(p[facei], TTR[facei]);
    }

    return tCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpvTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpvTR = volScalarField::New
    (
        "CpvTR",
        IOobject::NO_REGISTER,
        mesh,
        dimEnergy/dimMass/dimTemperature
    );
    auto& CpvTR = tCpvTR.ref();

    forAll(this->TTR_, celli)
    {
        CpvTR[celli] =
            this->cellMixture(celli).CpvTR(this->p_[celli], this->TTR_[celli]);
    }

    volScalarField::Boundary& CpvTRBf = CpvTR.boundaryFieldRef();

    forAll(CpvTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpvTR = CpvTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpvTR[facei] =
                this->patchFaceMixture(patchi, facei).CpvTR(pp[facei], pTTR[facei]);
        }
    }

    return tCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::CpByCpvTR
(
    const scalarField& p,
    const scalarField& TTR,
    const label patchi
) const
{
    auto tCpByCpvTR = tmp<scalarField>::New(TTR.size());
    auto& CpByCpvTR = tCpByCpvTR.ref();

    forAll(TTR, facei)
    {
        CpByCpvTR[facei] =
            this->patchFaceMixture(patchi, facei).CpByCpvTR(p[facei], TTR[facei]);
    }

    return tCpByCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::CpByCpvTR() const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tCpByCpvTR = volScalarField::New
    (
        "CpByCpvTR",
        IOobject::NO_REGISTER,
        mesh,
        dimless
    );
    auto& CpByCpvTR = tCpByCpvTR.ref();

    forAll(this->TTR_, celli)
    {
        CpByCpvTR[celli] = this->cellMixture(celli).CpByCpvTR
        (
            this->p_[celli],
            this->TTR_[celli]
        );
    }

    volScalarField::Boundary& CpByCpvTRBf =
        CpByCpvTR.boundaryFieldRef();

    forAll(CpByCpvTRBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pTTR = this->TTR_.boundaryField()[patchi];
        fvPatchScalarField& pCpByCpvTR = CpByCpvTRBf[patchi];

        forAll(pTTR, facei)
        {
            pCpByCpvTR[facei] = this->patchFaceMixture(patchi, facei).CpByCpvTR
            (
                pp[facei],
                pTTR[facei]
            );
        }
    }

    return tCpByCpvTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::THE_TR
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const labelList& cells
) const
{
    auto tTTR = tmp<scalarField>::New(h.size());
    auto& TTR = tTTR.ref();

    forAll(h, celli)
    {
        TTR[celli] =
            this->cellMixture(cells[celli]).THE_TR(h[celli], p[celli], TTR0[celli]);
    }

    return tTTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::THE_TR
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TTR0,
    const label patchi
) const
{

    auto tTTR = tmp<scalarField>::New(h.size());
    auto& TTR = tTTR.ref();

    forAll(h, facei)
    {
        TTR[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).THE_TR(h[facei], p[facei], TTR0[facei]);
    }

    return tTTR;
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TEH_Vib
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TVib0,
    const labelList& cells
) const
{

    auto tTVib = tmp<scalarField>::New(h.size());
    auto& TVib = tTVib.ref();

    forAll(h, celli)
    {
        TVib[celli] =
            this->cellMixture(cells[celli]).TEH_Vib(h[celli], p[celli], TVib0[celli]);
    }

    return tTVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::TEH_Vib
(
    const scalarField& h,
    const scalarField& p,
    const scalarField& TVib0,
    const label patchi
) const
{

    auto tTVib = tmp<scalarField>::New(h.size());
    auto& TVib = tTVib.ref();

    forAll(h, facei)
    {
        TVib[facei] = this->patchFaceMixture
        (
            patchi,
            facei
        ).TEH_Vib(h[facei], p[facei], TVib0[facei]);
    }

    return tTVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::W_s
(
) const
{
    const fvMesh& mesh = this->TTR_.mesh();

    auto tW_s = volScalarField::New
    (
        "W_s",
        IOobject::NO_REGISTER,
        mesh,
        dimMass/dimMoles
    );
    auto& W_s = tW_s.ref();

    scalarField& W_sCells = W_s.primitiveFieldRef();

    forAll(W_sCells, celli)
    {
        W_sCells[celli] = this->cellMixture(celli).W_s();
    }

    auto& W_sBf = W_s.boundaryFieldRef();

    forAll(W_sBf, patchi)
    {
        scalarField& W_sp = W_sBf[patchi];
        forAll(W_sp, facei)
        {
            W_sp[facei] = this->patchFaceMixture(patchi, facei).W_s();
        }
    }

    return tW_s;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaTR() const
{
    tmp<Foam::volScalarField> kappaTR(CpTR()*this->alpha_);
    kappaTR.ref().rename("kappaTR");
    return kappaTR;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaTR
(
    const label patchi
) const
{
    return
        CpTR
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}

template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaVib() const
{
    tmp<Foam::volScalarField> kappaVib(CpTR()*this->alpha_);
    kappaVib.ref().rename("kappaVib");
    return kappaVib;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField> Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaVib
(
    const label patchi
) const
{
    return
        CpTR
        (
            this->p_.boundaryField()[patchi],
            this->TVib_.boundaryField()[patchi],
            patchi
        )*this->alpha_.boundaryField()[patchi];
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaheTR() const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpvTR()*this->alpha_);
    alphaEff.ref().rename("alphaheTR");
    return alphaEff;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaheTR(const label patchi) const
{
    return
    this->CpByCpvTR
    (
        this->p_.boundaryField()[patchi],
        this->TTR_.boundaryField()[patchi],
        patchi
    )
   *this->alpha_.boundaryField()[patchi];
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> kappaEff(CpTR()*(this->alpha_ + alphat));
    kappaEff.ref().rename("kappaEff");
    return kappaEff;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::kappaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
        CpTR
        (
            this->p_.boundaryField()[patchi],
            this->TTR_.boundaryField()[patchi],
            patchi
        )
       *(
           this->alpha_.boundaryField()[patchi]
         + alphat
        );
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::volScalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaEff
(
    const volScalarField& alphat
) const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpvTR()*(this->alpha_ + alphat));
    alphaEff.ref().rename("alphaEff");
    return alphaEff;
}


template<class Basic2TThermo, class MixtureType>
Foam::tmp<Foam::scalarField>
Foam::he2TThermo<Basic2TThermo, MixtureType>::alphaEff
(
    const scalarField& alphat,
    const label patchi
) const
{
    return
    this->CpByCpvTR
    (
        this->p_.boundaryField()[patchi],
        this->TTR_.boundaryField()[patchi],
        patchi
    )
   *(
        this->alpha_.boundaryField()[patchi]
      + alphat
    );
}


template<class Basic2TThermo, class MixtureType>
bool Foam::he2TThermo<Basic2TThermo, MixtureType>::read()
{
    if (Basic2TThermo::read())
    {
        MixtureType::read(*this);
        return true;
    }

    return false;
}


// ************************************************************************* //
