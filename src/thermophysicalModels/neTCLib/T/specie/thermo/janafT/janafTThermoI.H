/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

#include "janafTThermo.H"
#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafTThermo<EquationOfState>::janafTThermo
(
    const EquationOfState& st,
    const scalar Tlow,
    const scalar Thigh,
    const scalar Tcommon,
    const typename janafTThermo<EquationOfState>::coeffArray& highCpCoeffs,
    const typename janafTThermo<EquationOfState>::coeffArray& lowCpCoeffs,
    const bool convertCoeffs
)
:
    EquationOfState(st),
    Tlow_(Tlow),
    Thigh_(Thigh),
    Tcommon_(Tcommon)
{
    if (convertCoeffs)
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel]*this->R();
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel]*this->R();
        }
    }
    else
    {
        for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
        {
            highCpCoeffs_[coefLabel] = highCpCoeffs[coefLabel];
            lowCpCoeffs_[coefLabel] = lowCpCoeffs[coefLabel];
        }
    }
}


template<class EquationOfState>
inline const typename Foam::janafTThermo<EquationOfState>::coeffArray&
Foam::janafTThermo<EquationOfState>::coeffs
(
    const scalar T
) const
{
    if (T < Tcommon_)
    {
        return lowCpCoeffs_;
    }
    else
    {
        return highCpCoeffs_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafTThermo<EquationOfState>::janafTThermo
(
    const word& name,
    const janafTThermo& jt
)
:
    EquationOfState(name, jt),
    Tlow_(jt.Tlow_),
    Thigh_(jt.Thigh_),
    Tcommon_(jt.Tcommon_)
{
    for (label coefLabel=0; coefLabel<nCoeffs_; coefLabel++)
    {
        highCpCoeffs_[coefLabel] = jt.highCpCoeffs_[coefLabel];
        lowCpCoeffs_[coefLabel] = jt.lowCpCoeffs_[coefLabel];
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    if (T < Tlow_ || T > Thigh_)
    {
        WarningInFunction
            << "attempt to use janafTThermo<EquationOfState>"
               " out of temperature range "
            << Tlow_ << " -> " << Thigh_ << ";  T = " << T
            << endl;

        return clamp(T, Tlow_, Thigh_);
    }

    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Tlow() const
{
    return Tlow_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Thigh() const
{
    return Thigh_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Tcommon() const
{
    return Tcommon_;
}


template<class EquationOfState>
inline const typename Foam::janafTThermo<EquationOfState>::coeffArray&
Foam::janafTThermo<EquationOfState>::highCpCoeffs() const
{
    return highCpCoeffs_;
}


template<class EquationOfState>
inline const typename Foam::janafTThermo<EquationOfState>::coeffArray&
Foam::janafTThermo<EquationOfState>::lowCpCoeffs() const
{
    return lowCpCoeffs_;
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        ((((a[4]*T + a[3])*T + a[2])*T + a[1])*T + a[0])
      + EquationOfState::Cp(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Ha
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        ((((a[4]/5.0*T + a[3]/4.0)*T + a[2]/3.0)*T + a[1]/2.0)*T + a[0])*T
      + a[5]
    ) + EquationOfState::H(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Hs
(
    const scalar p,
    const scalar T
) const
{
    return Ha(p, T) - Hc();
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Hc() const
{
    const coeffArray& a = lowCpCoeffs_;
    return
    (
        (
            (((a[4]/5.0*Tstd + a[3]/4.0)*Tstd + a[2]/3.0)*Tstd + a[1]/2.0)*Tstd
          + a[0]
        )*Tstd + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::S
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (((a[4]/4.0*T + a[3]/3.0)*T + a[2]/2.0)*T + a[1])*T + a[0]*log(T)
      + a[6]
    ) + EquationOfState::S(p, T);
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::Gstd
(
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
    (
        (
            a[0]*(1 - log(T))
          - (((a[4]/20.0*T + a[3]/12.0)*T + a[2]/6.0)*T + a[1]/2.0)*T
          - a[6]
        )*T
      + a[5]
    );
}


template<class EquationOfState>
inline Foam::scalar Foam::janafTThermo<EquationOfState>::dCpdT
(
    const scalar p,
    const scalar T
) const
{
    const coeffArray& a = coeffs(T);
    return
        (((4*a[4]*T + 3*a[3])*T + 2*a[2])*T + a[1]);
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline void Foam::janafTThermo<EquationOfState>::operator+=
(
    const janafTThermo<EquationOfState>& jt
)
{
    scalar Y1 = this->Y();

    EquationOfState::operator+=(jt);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = jt.Y()/this->Y();

        Tlow_ = max(Tlow_, jt.Tlow_);
        Thigh_ = min(Thigh_, jt.Thigh_);

        if
        (
            janafTThermo<EquationOfState>::debug
         && !equal(Tcommon_, jt.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << Tcommon_ << " for "
                << (this->name().size() ? this->name() : "others")
                << " != " << jt.Tcommon_ << " for "
                << (jt.name().size() ? jt.name() : "others")
                << exit(FatalError);
        }

        for
        (
            label coefLabel=0;
            coefLabel<janafTThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs_[coefLabel] =
                Y1*highCpCoeffs_[coefLabel]
              + Y2*jt.highCpCoeffs_[coefLabel];

            lowCpCoeffs_[coefLabel] =
                Y1*lowCpCoeffs_[coefLabel]
              + Y2*jt.lowCpCoeffs_[coefLabel];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::janafTThermo<EquationOfState> Foam::operator+
(
    const janafTThermo<EquationOfState>& jt1,
    const janafTThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs = jt1;
    eofs += jt2;

    if (mag(eofs.Y()) < SMALL)
    {
        return janafTThermo<EquationOfState>
        (
            eofs,
            jt1.Tlow_,
            jt1.Thigh_,
            jt1.Tcommon_,
            jt1.highCpCoeffs_,
            jt1.lowCpCoeffs_
        );
    }
    else
    {
        const scalar Y1 = jt1.Y()/eofs.Y();
        const scalar Y2 = jt2.Y()/eofs.Y();

        typename janafTThermo<EquationOfState>::coeffArray highCpCoeffs;
        typename janafTThermo<EquationOfState>::coeffArray lowCpCoeffs;

        for
        (
            label coefLabel=0;
            coefLabel<janafTThermo<EquationOfState>::nCoeffs_;
            coefLabel++
        )
        {
            highCpCoeffs[coefLabel] =
                Y1*jt1.highCpCoeffs_[coefLabel]
              + Y2*jt2.highCpCoeffs_[coefLabel];

            lowCpCoeffs[coefLabel] =
                Y1*jt1.lowCpCoeffs_[coefLabel]
              + Y2*jt2.lowCpCoeffs_[coefLabel];
        }

        if
        (
            janafTThermo<EquationOfState>::debug
         && !equal(jt1.Tcommon_, jt2.Tcommon_)
        )
        {
            FatalErrorInFunction
                << "Tcommon " << jt1.Tcommon_ << " for "
                << (jt1.name().size() ? jt1.name() : "others")
                << " != " << jt2.Tcommon_ << " for "
                << (jt2.name().size() ? jt2.name() : "others")
                << exit(FatalError);
        }

        return janafTThermo<EquationOfState>
        (
            eofs,
            max(jt1.Tlow_, jt2.Tlow_),
            min(jt1.Thigh_, jt2.Thigh_),
            jt1.Tcommon_,
            highCpCoeffs,
            lowCpCoeffs
        );
    }
}


template<class EquationOfState>
inline Foam::janafTThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const janafTThermo<EquationOfState>& jt
)
{
    return janafTThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(jt),
        jt.Tlow_,
        jt.Thigh_,
        jt.Tcommon_,
        jt.highCpCoeffs_,
        jt.lowCpCoeffs_
    );
}


template<class EquationOfState>
inline Foam::janafTThermo<EquationOfState> Foam::operator==
(
    const janafTThermo<EquationOfState>& jt1,
    const janafTThermo<EquationOfState>& jt2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(jt1)
     == static_cast<const EquationOfState&>(jt2)
    );

    const scalar Y1 = jt2.Y()/eofs.Y();
    const scalar Y2 = jt1.Y()/eofs.Y();

    typename janafTThermo<EquationOfState>::coeffArray highCpCoeffs;
    typename janafTThermo<EquationOfState>::coeffArray lowCpCoeffs;

    for
    (
        label coefLabel=0;
        coefLabel<janafTThermo<EquationOfState>::nCoeffs_;
        coefLabel++
    )
    {
        highCpCoeffs[coefLabel] =
            Y1*jt2.highCpCoeffs_[coefLabel]
          - Y2*jt1.highCpCoeffs_[coefLabel];

        lowCpCoeffs[coefLabel] =
            Y1*jt2.lowCpCoeffs_[coefLabel]
          - Y2*jt1.lowCpCoeffs_[coefLabel];
    }

    if
    (
        janafTThermo<EquationOfState>::debug
     && !equal(jt2.Tcommon_, jt1.Tcommon_)
    )
    {
        FatalErrorInFunction
            << "Tcommon " << jt2.Tcommon_ << " for "
            << (jt2.name().size() ? jt2.name() : "others")
            << " != " << jt1.Tcommon_ << " for "
            << (jt1.name().size() ? jt1.name() : "others")
            << exit(FatalError);
    }

    return janafTThermo<EquationOfState>
    (
        eofs,
        max(jt2.Tlow_, jt1.Tlow_),
        min(jt2.Thigh_, jt1.Thigh_),
        jt2.Tcommon_,
        highCpCoeffs,
        lowCpCoeffs
    );
}


// ************************************************************************* //
