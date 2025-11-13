/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "WilkeMR.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType>
Foam::WilkeMR<MixtureType>::WilkeMR
(
    const MixtureType& mixture
)
:
    baseMR<MixtureType>(mixture),
    mix_(mixture),
    names_(mixture.species()),
    Xi_
    (
        IOobject
        (
            "Xi",
            mixture.Y()[0].time().timeName(),  // instance: current time
            mixture.Y()[0].mesh(),             // registry: mesh (fvMesh is an objectRegistry)
            Foam::IOobject::NO_READ,
            Foam::IOobject::NO_WRITE,
            false                           // do not register in DB (scratch field)
        ),
        mixture.Y()[0].mesh(),
        Foam::dimensionedScalar("zero", Foam::dimless, 0.0)
    )
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
template<class MixtureType>
Foam::volScalarField&
Foam::WilkeMR<MixtureType>::computeXiFromYi
(
    const label speciei
)
{
    const PtrList<volScalarField>& Y = mix_.Y();
    const volScalarField& Yi = Y[speciei];
    const scalar Wi = mix_.W(speciei);

    // Denominator: sum_j ( Yj / Wj )
    // Start from a small positive value to avoid divide-by-zero.
    tmp<volScalarField> tDenom
    (
        new volScalarField
        (
            IOobject
            (
                "Xi_denom",
                Yi.time().timeName(),
                Yi.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            Yi.mesh(),
            dimensionedScalar("small", dimless, SMALL)
        )
    );

    volScalarField& denom = tDenom.ref();

    for (label j = 0; j < Y.size(); ++j)
    {
        const scalar Wj = mix_.W(j);
        denom += Y[j]/Wj;            // denom += Y_j / W_j
    }

    // Xi_i = (Yi / Wi) / denom
    Xi_ = (Yi/Wi)/denom;

    return Xi_;
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::scalingFactor
(
    const label speciei,      // specie index speciei
    const label celli,  // cell index
    const scalar TTR      // local temperature
)
{
    const label nSpec = mix_.Y().size();

    // --- X_s from existing helper ---
    volScalarField& Xi_s = computeXiFromYi(speciei);
    scalar Xs = Xi_s[celli];

    const scalar Ms   = mix_.W(speciei);
    const scalar mu_s = mix_.mu(speciei, TTR);

    // start φ_s with X_s
    scalar phi_s = Xs;

    // --- sum over r ≠ speciei ---
    for (label r = 0; r < nSpec; ++r)
    {
        if (r == speciei) continue;

        // X_r from same helper (Xi_ is reused internally)
        volScalarField& Xi_r = computeXiFromYi(r);
        scalar Xr = Xi_r[celli];

        const scalar Mr   = mix_.W(r);
        const scalar mu_r = mix_.mu(r, TTR);

        // [1 + sqrt(mu_s/mu_r) * (Mr/Ms)^(1/4)]^2
        const scalar term1 =
            1.0 + Foam::sqrt(mu_s/mu_r)*Foam::pow(Mr/Ms, 0.25);
        const scalar term1sq = Foam::sqr(term1);

        // [ sqrt( 8 (1 + Ms/Mr) ) ]^-1
        const scalar term2 =
            1.0/Foam::sqrt(8.0*(1.0 + Ms/Mr));

        phi_s += Xr * term1sq * term2;
    }

    return phi_s;
}

template<class MixtureType>
template<class QGetter>
Foam::scalar
Foam::WilkeMR<MixtureType>::QCell
(
    const label celli,
    const scalar p,
    const scalar TTR,
    const QGetter& getQ
)
{
    const PtrList<volScalarField>& Y = mix_.Y();
    const label nSpec = Y.size();

    scalar Qmix = 0.0;          

    for (label speciei = 0; speciei < nSpec; ++speciei)
    {

        volScalarField& Xi_s = this->computeXiFromYi(speciei);
        const scalar Xs = Xi_s[celli];

        const scalar phi_s = this->scalingFactor(speciei, celli, TTR);

        const scalar Qs = getQ(speciei, p, TTR);

        Qmix += Qs*Xs/phi_s;
    }

    return Qmix;
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::muCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getMu = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.mu(speciei, TTRi);       
    };

    return QCell(celli, p, TTR, getMu);
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::kappaTRCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getK = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.kappaTR(speciei, p, TTRi);    
    };

    return QCell(celli, p, TTR, getK);
}

template<class MixtureType>
Foam::scalar
Foam::WilkeMR<MixtureType>::alphaTRCell
(
    const label celli,
    const scalar p,
    const scalar TTR
)
{
    auto getAlpha = [&](const label speciei, const scalar p, const scalar TTRi)
    {
        return mix_.kappaTR(speciei, p, TTRi)/mix_.CpTR(speciei, p, TTRi);    
    };

    return QCell(celli, p, TTR, getAlpha);
}

// ************************************************************************* //
