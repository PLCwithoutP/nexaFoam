/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2017 OpenFOAM Foundation
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

Class
    Foam::sensibleEnthalpyT

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleEnthalpyT_H
#define sensibleEnthalpyT_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensibleEnthalpyT Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoT>
class sensibleEnthalpyT
{

public:

    // Constructors

        //- Construct
        sensibleEnthalpyT()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleEnthalpyT";
        }

        // Fundamental properties

            static word energyName()
            {
                return "h";
            }

            //- Heat capacity at constant pressure [J/(kg K)]
            scalar Cpv
            (
                const ThermoT& thermoT,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar cp = thermoT.Cp(p, T);
                return cp;
                #else
                return thermoT.Cp(p, T);
                #endif
            }

            //- Cp/Cp []
            scalar CpByCpv
            (
                const ThermoT& thermoT,
                const scalar p,
                const scalar T
            ) const
            {
                return 1;
            }

            //- Sensible enthalpy [J/kg]
            scalar HE
            (
                const ThermoT& thermoT,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar hs = thermoT.Hs(p, T);
                return hs;
                #else
                return thermoT.Hs(p, T);
                #endif
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar THE
            (
                const ThermoT& thermoT,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar ths = thermoT.THs(h, p, T0);
                return ths;
                #else
                return thermoT.THs(h, p, T0);
                #endif
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
