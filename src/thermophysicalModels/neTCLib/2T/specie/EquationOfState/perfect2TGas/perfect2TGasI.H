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

#include "perfect2TGas.H"
#include "thermodynamicConstants.H"

using namespace Foam::constant::thermodynamic;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2TGas<Specie>::perfect2TGas(const Specie& sp)
:
    Specie(sp)
{}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2TGas<Specie>::perfect2TGas
(
    const word& name,
    const perfect2TGas<Specie>& pg
)
:
    Specie(name, pg)
{}


template<class Specie>
inline Foam::autoPtr<Foam::perfect2TGas<Specie>>
Foam::perfect2TGas<Specie>::clone() const
{
    return autoPtr<perfect2TGas<Specie>>::New(*this);
}


template<class Specie>
inline Foam::autoPtr<Foam::perfect2TGas<Specie>> Foam::perfect2TGas<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<perfect2TGas<Specie>>::New(dict);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::rho(scalar p, scalar T) const
{
    return p/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::HT(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::CpT(scalar p, scalar T) const
{
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::ET(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::ER(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::EV(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::CvT(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::CvR(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::CvVib(scalar p, scalar T) const
{
    return 0;
}

template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::S(scalar p, scalar T) const
{
    return -this->R()*log(p/Pstd);
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::psi(scalar p, scalar T) const
{
    return 1.0/(this->R()*T);
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::Z(scalar p, scalar T) const
{
    return 1;
}


template<class Specie>
inline Foam::scalar Foam::perfect2TGas<Specie>::CpMCv(scalar p, scalar T) const
{
    return this->R();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Specie>
inline void Foam::perfect2TGas<Specie>::operator+=(const perfect2TGas<Specie>& pg)
{
    Specie::operator+=(pg);
}


template<class Specie>
inline void Foam::perfect2TGas<Specie>::operator*=(const scalar s)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::perfect2TGas<Specie> Foam::operator+
(
    const perfect2TGas<Specie>& pg1,
    const perfect2TGas<Specie>& pg2
)
{
    return perfect2TGas<Specie>
    (
        static_cast<const Specie&>(pg1) + static_cast<const Specie&>(pg2)
    );
}


template<class Specie>
inline Foam::perfect2TGas<Specie> Foam::operator*
(
    const scalar s,
    const perfect2TGas<Specie>& pg
)
{
    return perfect2TGas<Specie>(s*static_cast<const Specie&>(pg));
}


template<class Specie>
inline Foam::perfect2TGas<Specie> Foam::operator==
(
    const perfect2TGas<Specie>& pg1,
    const perfect2TGas<Specie>& pg2
)
{
    return perfect2TGas<Specie>
    (
        static_cast<const Specie&>(pg1) == static_cast<const Specie&>(pg2)
    );
}


// ************************************************************************* //
