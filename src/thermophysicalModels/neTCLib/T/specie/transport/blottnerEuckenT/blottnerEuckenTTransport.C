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

#include "blottnerEuckenTTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoT>
Foam::scalar Foam::blottnerEuckenTTransport<ThermoT>::readCoeff
(
    const word& coeffName,
    const dictionary& dict
)
{
    return dict.subDict("transport").get<scalar>(coeffName);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoT>
Foam::blottnerEuckenTTransport<ThermoT>::blottnerEuckenTTransport(const dictionary& dict)
:
    ThermoT(dict),
    AB_(readCoeff("AB", dict)),
    BB_(readCoeff("BB", dict)),
    CB_(readCoeff("CB", dict))
{}


template<class ThermoT>
Foam::blottnerEuckenTTransport<ThermoT>::blottnerEuckenTTransport
(
    const ThermoT& t,
    const dictionary& dict
)
:
    ThermoT(t),
    AB_(readCoeff("AB", dict)),
    BB_(readCoeff("BB", dict)),
    CB_(readCoeff("CB", dict))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoT>
void Foam::blottnerEuckenTTransport<ThermoT>::write(Ostream& os) const
{
    os.beginBlock(this->specie::name());

    ThermoT::write(os);

    // Entries in dictionary format
    {
        os.beginBlock("transport");
        os.writeEntry("AB", AB_);
        os.writeEntry("BB", BB_);
        os.writeEntry("CB", CB_);
        os.endBlock();
    }

    os.endBlock();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class ThermoT>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const blottnerEuckenTTransport<ThermoT>& st
)
{
    st.write(os);
    return os;
}


// ************************************************************************* //
