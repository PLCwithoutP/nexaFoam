/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
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
    Foam::typeDecider

Group
    grp

Description
    Decide which temperature model (1T or 2T) to use.

    Template params:
    - neTThermo : the 1-temperature thermo class (e.g. neTThermo)
    - ne2TThermo : the 2-temperature thermo class (e.g. ne2TThermo)
SourceFiles
    typeDeciderI.H
    typeDecider.C

\*---------------------------------------------------------------------------*/

#ifndef typeDecider_H
#define typeDecider_H

#include "autoPtr.H"
#include "word.H"
#include "fvMesh.H"
#include "error.H"

namespace Foam
{

template<class OneT, class TwoT>
class typeDecider
{
    protected:
    
        IOdictionary thermoDict;
        const dictionary tt;
        const word tempModel;

    public:

    // Constructor
        //- Construct from mesh and phase name
        typeDecider
        (
            const fvMesh&
        );

    // Member Functions
    
        // - Construct exactly one of the thermos based on `type`
        bool construct
        (
            autoPtr<OneT>& p1T,
            autoPtr<TwoT>& p2T,
            const fvMesh& mesh
        );
};

} // namespace Foam

#include "typeDeciderI.H"

#ifdef NoRepository
#   include "typeDecider.C"
#endif

#endif
