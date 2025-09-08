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

#include "neTThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janafTThermo.H"
#include "sensibleEnthalpyT.H"
#include "sensibleInternalEnergyT.H"
#include "thermoT.H"

#include "heNeTThermo.H"
#include "pureMixture.H"

#include "thermoPhysicsTTypes.H"

#include "blottnerEuckenTTransport.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeThermos
(
    neTThermo,
    heNeTThermo,
    pureMixture,
    blottnerEuckenTTransport,
    sensibleEnthalpyT,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    neTThermo,
    heNeTThermo,
    pureMixture,
    blottnerEuckenTTransport,
    sensibleEnthalpyT,
    janafTThermo,
    perfectGas,
    specie
);



/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */


makeThermos
(
    neTThermo,
    heNeTThermo,
    pureMixture,
    blottnerEuckenTTransport,
    sensibleInternalEnergyT,
    eConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    neTThermo,
    heNeTThermo,
    pureMixture,
    blottnerEuckenTTransport,
    sensibleInternalEnergyT,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    neTThermo,
    heNeTThermo,
    pureMixture,
    blottnerEuckenTTransport,
    sensibleInternalEnergyT,
    janafTThermo,
    perfectGas,
    specie
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
