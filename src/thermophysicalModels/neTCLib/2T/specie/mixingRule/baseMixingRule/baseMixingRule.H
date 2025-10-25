/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014 OpenFOAM Foundation
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
    Foam::baseMixingRule

Description
    Base-class for all mixing rules in order to calculate mixture quantities.

SourceFiles
    baseMixingRule.C

\*---------------------------------------------------------------------------*/

#ifndef baseMixingRule_H
#define baseMixingRule_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration
class basicMultiComponent2TMixture;

/*---------------------------------------------------------------------------*\
                Class baseMixingRule Declaration
\*---------------------------------------------------------------------------*/

class baseMixingRule
{
private:
    // Private Member Functions

        //- No copy construct
        baseMixingRule(const baseMixingRule&) = delete;

        //- No copy assignment
        void operator=(const baseMixingRule&) = delete;

protected:

    // Protected member variables

    // Must be multiComponent2TMixture, otherwise will fail
    const basicMultiComponent2TMixture& mixture_;

public:

    //- Runtime type information
    TypeName("baseMixingRule");


    // Constructors

        //- Construct from components
        baseMixingRule();


    //- Destructor
    virtual ~baseMixingRule();


    // Member Functions

        //- Return the scalar quantity
        virtual tmp<volScalarField> Qscalar(const volScalarField& Yi) const = 0;

        //- Return the vectoral quantity
        virtual tmp<volVectorField> Qvector(const volScalarField& Yi) const = 0;

        //- Return the scalar quantity for patch
        virtual tmp<scalarField> Qscalar(const label patchi, const volScalarField& Yi) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //