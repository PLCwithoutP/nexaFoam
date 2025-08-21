/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "specie.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::blottnerEuckenTransport<Thermo>::blottnerEuckenTransport
(
    const Thermo& t,
    const scalar AB,
    const scalar BB,
    const scalar CB
)
:
    Thermo(t),
    AB_(AB),
    BB_(BB),
    CB_(CB)
{}


template<class Thermo>
inline Foam::blottnerEuckenTransport<Thermo>::blottnerEuckenTransport
(
    const word& name,
    const blottnerEuckenTransport& st
)
:
    Thermo(name, st),
    AB_(st.AB_),
    BB_(st.BB_),
    CB_(st.CB_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::blottnerEuckenTransport<Thermo>>
Foam::blottnerEuckenTransport<Thermo>::clone() const
{
    return autoPtr<blottnerEuckenTransport<Thermo>>::New(*this);
}


template<class Thermo>
inline Foam::autoPtr<Foam::blottnerEuckenTransport<Thermo>>
Foam::blottnerEuckenTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<blottnerEuckenTransport<Thermo>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::blottnerEuckenTransport<Thermo>::mu
(
    const scalar T
) const
{
    //Info << "Function called" << endl;
    return 0.1*exp((AB_*log(T) + BB_)*log(T) + CB_);
}


template<class Thermo>
inline Foam::scalar Foam::blottnerEuckenTransport<Thermo>::kappa_TR
(
    const scalar p, const scalar T
) const
{
    Info << "kappa_TR function called" << endl;
    scalar Cv_T_ = this->Cv(p, T);
    scalar Cv_R_ = this->Cv(p, T);
    return mu(T)*((5/2)*Cv_T_ + Cv_R_);
}

template<class Thermo>
inline Foam::scalar Foam::blottnerEuckenTransport<Thermo>::kappa_V
(
    const scalar T
) const
{
    scalar nv_ = 1.2;
    return nv_*mu(T);
}


template<class Thermo>
inline Foam::scalar Foam::blottnerEuckenTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{

    return kappa_TR(p, T)/this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::blottnerEuckenTransport<Thermo>::operator+=
(
    const blottnerEuckenTransport<Thermo>& st
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(st);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = st.Y()/this->Y();

        AB_ = Y1*AB_ + Y2*st.AB_;
        BB_ = Y1*BB_ + Y2*st.BB_;
        CB_ = Y1*CB_ + Y2*st.CB_;
    }
}


template<class Thermo>
inline void Foam::blottnerEuckenTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::blottnerEuckenTransport<Thermo> Foam::operator+
(
    const blottnerEuckenTransport<Thermo>& st1,
    const blottnerEuckenTransport<Thermo>& st2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(st1) + static_cast<const Thermo&>(st2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return blottnerEuckenTransport<Thermo>
        (
            t,
            st1.AB_,
            st1.BB_,
            st1.CB_
        );
    }
    else
    {
        scalar Y1 = st1.Y()/t.Y();
        scalar Y2 = st2.Y()/t.Y();

        return blottnerEuckenTransport<Thermo>
        (
            t,
            Y1*st1.AB_ + Y2*st2.AB_,
            Y1*st1.BB_ + Y2*st2.BB_,
            Y1*st1.CB_ + Y2*st2.CB_
        );
    }
}


template<class Thermo>
inline Foam::blottnerEuckenTransport<Thermo> Foam::operator*
(
    const scalar s,
    const blottnerEuckenTransport<Thermo>& st
)
{
    return blottnerEuckenTransport<Thermo>
    (
        s*static_cast<const Thermo&>(st),
        st.AB_,
        st.BB_,
        st.CB_
    );
}


// ************************************************************************* //
