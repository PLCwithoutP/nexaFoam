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

#include "VTEnergySource.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::VTEnergySource<MixtureType, MixingRule>::VTEnergySource
(
    const MixtureType& mixture,
    const MixingRule& MR
)
:
    baseEnergySource<MixtureType, MixingRule>(mixture),
    mix_(mixture),
    mr_(MR),
    names_(mixture.species()),
    sigma_prime_(3e-21), // N2, O2, NO
    boltzmann_const_(1.380649e-23)
{

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::sigmai
(
    const scalar TTR 
)
{
    return sigma_prime_*(50000.0/TTR)*(50000.0/TTR);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::c_bar_s
(
    const scalar TTR,
    const label s 
)
{
    return pow((8*mix_.R()*TTR)/(constant::mathematical::pi*Wi(s)) , 0.5);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::n_s
(
    const scalar p,
    const label TTR
)
{
    return p/(TTR*boltzmann_const_);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_P_sr
(
    const scalar p,
    const scalar TTR,
    const label s
)
{
    return 1/(c_bar_s(TTR,s) * sigmai(TTR) * n_s(p, TTR));
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::A_sr
(
    const label s,
    const label r
)
{
    scalar val = 1.16 * pow(10.0, -9.0/2.0);
    return val * pow((Wi(s)*Wi(r))/(Wi(s) + Wi(r)) , 0.5) * pow(thetai(s), 4.0/3.0);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::B_sr
(
    const label s,
    const label r
)
{
    scalar val = 0.015 * pow(10.0, -3.0/4.0);
    return val * pow((Wi(s)*Wi(r))/(Wi(s) + Wi(r)) , 0.25);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_MW_sr
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label r
)
{
    scalar A = A_sr(s,r);
    scalar B = B_sr(s,r);
    return (1/p)*exp(A * (pow(TTR, -1/3) - B) - 18.42);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_sr
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label r
)
{
    return T_MW_sr(p, TTR, s, r) + T_P_sr(p, TTR, s);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::T_s
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label celli
)
{
    scalar sumX          = 0.0;
    scalar sumX_over_tau = 0.0;

    for (label r = 0; r < nSpecie_; ++r)
    {
        if (!isMolecular(r)) continue;  // corresponds to r = mol in the sum

        // X_r at this cell
        const volScalarField& XrField = mr_.computeXiFromYi(r);
        const scalar Xr = XrField[celli];

        if (Xr <= SMALL) continue;     // skip negligible species

        const scalar tau_sr = T_sr(p, TTR, s, r);

        sumX          += Xr;
        sumX_over_tau += Xr/tau_sr;
    }

    // safety: no valid colliders → avoid 0/0
    if (sumX_over_tau <= SMALL)
    {
        // choose what makes sense for you:
        // return GREAT; or 0.0; or FatalErrorInFunction << ... << abort(FatalError);
        return GREAT;
    }

    // τ_{s,V-T} = Σ X_r / Σ (X_r / τ_{s-r})
    return sumX/sumX_over_tau;
}

// Classical Landau-Teller Formula
template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VTEnergySource<MixtureType, MixingRule>::Q_TR_s
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const label s,
    const label celli
)
{
    return mix_.rho(s,p,TTR)*(mix_.EsVib(s,p,TTR,TTR,thetai(s)) - mix_.EsVib(s,p,TTR,TVib,thetai(s)))/(T_s(p, TTR, s, celli));
}
// ************************************************************************* //
