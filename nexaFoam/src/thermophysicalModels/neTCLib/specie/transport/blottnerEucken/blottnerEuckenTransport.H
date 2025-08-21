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
    Foam::blottnerEuckenTransport

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
    blottnerEuckenTransportI.H
    blottnerEuckenTransport.C

\*---------------------------------------------------------------------------*/

#ifndef blottnerEuckenTransport_H
#define blottnerEuckenTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations

template<class Thermo> class blottnerEuckenTransport;

template<class Thermo>
inline blottnerEuckenTransport<Thermo> operator+
(
    const blottnerEuckenTransport<Thermo>&,
    const blottnerEuckenTransport<Thermo>&
);

template<class Thermo>
inline blottnerEuckenTransport<Thermo> operator*
(
    const scalar,
    const blottnerEuckenTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const blottnerEuckenTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                    Class blottnerEuckenTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class blottnerEuckenTransport
:
    public Thermo
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
        inline blottnerEuckenTransport
        (
            const Thermo& t,
            const scalar AB,
            const scalar BB,
            const scalar CB
        );

        //- Construct as named copy
        inline blottnerEuckenTransport(const word&, const blottnerEuckenTransport&);

        //- Construct from dictionary
        explicit blottnerEuckenTransport(const dictionary& dict);

        //- Construct from base thermo and dictionary
        blottnerEuckenTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<blottnerEuckenTransport> clone() const;

        // Selector from dictionary
        inline static autoPtr<blottnerEuckenTransport> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "blottnerEucken<" + Thermo::typeName() + '>';
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

        inline void operator+=(const blottnerEuckenTransport&);

        inline void operator*=(const scalar);


    // Friend Operators

        friend blottnerEuckenTransport operator+ <Thermo>
        (
            const blottnerEuckenTransport&,
            const blottnerEuckenTransport&
        );

        friend blottnerEuckenTransport operator* <Thermo>
        (
            const scalar,
            const blottnerEuckenTransport&
        );


    // IOstream Operators

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const blottnerEuckenTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blottnerEuckenTransportI.H"

#ifdef NoRepository
    #include "blottnerEuckenTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
