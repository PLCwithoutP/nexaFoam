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

#include "thermoT.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class ThermoT, template<class> class Type>
const Foam::scalar Foam::species::thermoT<ThermoT, Type>::tol_ = 1.0e-4;

template<class ThermoT, template<class> class Type>
const int Foam::species::thermoT<ThermoT, Type>::maxIter_ = 100;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
Foam::species::thermoT<ThermoT, Type>::thermoT(const dictionary& dict)
:
    ThermoT(dict)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
void Foam::species::thermoT<ThermoT, Type>::write(Ostream& os) const
{
    ThermoT::write(os);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
Foam::Ostream& Foam::species::operator<<
(
    Ostream& os, const thermoT<ThermoT, Type>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
