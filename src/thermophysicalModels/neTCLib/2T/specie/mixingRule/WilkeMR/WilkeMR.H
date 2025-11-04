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
    Foam::WilkeMR

Group
    grpSpecieMixingRule

Description
    Wilke's mixing rule class using mixture class.

SourceFiles
    WilkeMRI.H
    WilkeMR.C

\*---------------------------------------------------------------------------*/

#ifndef WilkeMR_H
#define WilkeMR_H

#include "baseMixingRule.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class WilkeMR Declaration
\*---------------------------------------------------------------------------*/


class WilkeMR
:
    public baseMixingRule
{
    public:

    private:
        // Private Data
            List<scalar> W_;       // [kg/kmol] or [g/mol]

            PtrList<volScalarField> mu_;  // species viscosities
            PtrList<volScalarField> x_;   // mole fractions
    
    private:
        //- No copy construct
        WilkeMR(const WilkeMR&) = delete;

        //- No copy assignment
        void operator=(const WilkeMR&) = delete;

    public:
        // Constructors

        // Member Functions

            //- Return the instantiated type name
            static word typeName()
            {
                return "WilkeMR";
            }

            //- Return the scalar quantity
            inline tmp<volScalarField> Qscalar(const volScalarField& Yi) const;

            //- Return the scalar quantity for patch
            inline tmp<scalarField> Qscalar(const label patchi, const volScalarField& Yi) const;

            //- Return the scaling factor for mixture
            inline scalar phiSpecie(const volScalarField& Yi) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
    #include "WilkeMR.C"
#endif

#endif //WilkeMR_H