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
    Foam::sensibleEnthalpy2T

Group
    grpSpecieThermo

Description
    Thermodynamics mapping class to expose the sensible enthalpy functions.

\*---------------------------------------------------------------------------*/

#ifndef sensibleEnthalpy2T_H
#define sensibleEnthalpy2T_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class sensibleEnthalpy2T Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo2T>
class sensibleEnthalpy2T
{

public:

    // Constructors

        //- Construct
        sensibleEnthalpy2T()
        {}


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "sensibleEnthalpy2T";
        }

        // Fundamental properties

            static word energyName()
            {
                return "h";
            }

            //- Heat capacity at constant pressure [J/(kg K)]
            scalar Cpv
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar cp = thermo2T.Cp(p, T);
                return cp;
                #else
                return thermo2T.Cp(p, T);
                #endif
            }

            //- Cp/Cp []
            scalar CpByCpv
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar T
            ) const
            {
                return 1;
            }

            //- Sensible enthalpy [J/kg]
            scalar HE
            (
                const Thermo2T& thermo2T,
                const scalar p,
                const scalar T
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar hs = thermo2T.Hs(p, T);
                return hs;
                #else
                return thermo2T.Hs(p, T);
                #endif
            }

            //- Temperature from sensible enthalpy
            //  given an initial temperature T0
            scalar THE
            (
                const Thermo2T& thermo2T,
                const scalar h,
                const scalar p,
                const scalar T0
            ) const
            {
                #ifdef __clang__
                // Using volatile to prevent compiler optimisations leading to
                // a sigfpe
                volatile const scalar ths = thermo2T.THs(h, p, T0);
                return ths;
                #else
                return thermo2T.THs(h, p, T0);
                #endif
            }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
