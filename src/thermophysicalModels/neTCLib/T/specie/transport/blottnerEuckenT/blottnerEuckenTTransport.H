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

Class
    Foam::blottnerEuckenTTransport

Group
    grpSpecieTransport

Description
    Transport package using Blottner's viscosity formula and Eucken's
    thermal conductivity formula.

    Templated into a given thermodynamics package (needed for thermal
    conductivity).

    Dynamic viscosity [kg/m.s]
    \f[
        \mu_s = 0.1 * exp((A_B*ln(T) + B_B)*ln(T) + C_B)
    \f]
    T is translational temperature.

    Thermal conductivity [kg/m.s]
    \f[
        \kappa_tr = \mu_s (\frac{5}{2}*c_{v,t,s} + c_{v,r,s})
        \kappa_v = \nu_v * \mu_s
        \nu_v = 1.2
    \f]
    Holds until 6000 K.

SourceFiles
    blottnerEuckenTTransportI.H
    blottnerEuckenTTransport.C

\*---------------------------------------------------------------------------*/

#ifndef blottnerEuckenTTransport_H
#define blottnerEuckenTTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class ThermoT> class blottnerEuckenTTransport;

template<class ThermoT>
inline blottnerEuckenTTransport<ThermoT> operator+
(
    const blottnerEuckenTTransport<ThermoT>&,
    const blottnerEuckenTTransport<ThermoT>&
);

template<class ThermoT>
inline blottnerEuckenTTransport<ThermoT> operator*
(
    const scalar,
    const blottnerEuckenTTransport<ThermoT>&
);

template<class ThermoT>
Ostream& operator<<
(
    Ostream&,
    const blottnerEuckenTTransport<ThermoT>&
);


/*---------------------------------------------------------------------------*\
                    Class blottnerEuckenTTransport Declaration
\*---------------------------------------------------------------------------*/

template<class ThermoT>
class blottnerEuckenTTransport
:
    public ThermoT
{
    // Private Data

        // Blottner's coefficients
        scalar AB_, BB_, CB_;


    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Generated Methods: copy construct, copy assignment


    // Constructors

        //- Construct from components
        inline blottnerEuckenTTransport
        (
            const ThermoT& t,
            const scalar AB,
            const scalar BB,
            const scalar CB
        );

        //- Construct as named copy
        inline blottnerEuckenTTransport(const word&, const blottnerEuckenTTransport&);

        //- Construct from dictionary
        explicit blottnerEuckenTTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        blottnerEuckenTTransport(const ThermoT& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<blottnerEuckenTTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<blottnerEuckenTTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "blottnerEuckenT<" + ThermoT::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar T) const;

        //- Thermal translational-rotational conductivity [W/mK]
        inline scalar kappa_TR(const scalar p, const scalar T) const;

        //- Thermal vibrational conductivity [W/mK]
        inline scalar kappa_V(const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        // Species diffusivity
        //inline scalar D(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const blottnerEuckenTTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend blottnerEuckenTTransport operator+ <ThermoT>
        (
            const blottnerEuckenTTransport&,
            const blottnerEuckenTTransport&
        );

        friend blottnerEuckenTTransport operator* <ThermoT>
        (
            const scalar,
            const blottnerEuckenTTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <ThermoT>
        (
            Ostream&,
            const blottnerEuckenTTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blottnerEuckenTTransportI.H"

#ifdef NoRepository
    #include "blottnerEuckenTTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
