/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2021,2022 OpenCFD Ltd.
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

#include "ne2TThermo.H"
#include "makeThermo.H"

#include "specie.H"
#include "perfectGas.H"
#include "PengRobinsonGas.H"
#include "hConstThermo.H"
#include "eConstThermo.H"
#include "janaf2TThermo.H"
#include "sensibleEnthalpy2T.H"
#include "sensibleInternalEnergy2T.H"
#include "thermo2T.H"

#include "hPolynomialThermo.H"

#include "heNe2TThermo.H"
#include "pureZoneMixture.H"

#include "thermoPhysics2TTypes.H"

#include "blottnerEucken2TTransport.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * * Enthalpy-based * * * * * * * * * * * * * */

makeThermos
(
    ne2TThermo,
    heNe2TThermo,
    pureZoneMixture,
    blottnerEucken2TTransport,
    sensibleEnthalpy2T,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    ne2TThermo,
    heNe2TThermo,
    pureZoneMixture,
    blottnerEucken2TTransport,
    sensibleEnthalpy2T,
    janaf2TThermo,
    perfectGas,
    specie
);




/* * * * * * * * * * * * * * Internal-energy-based * * * * * * * * * * * * * */


makeThermos
(
    ne2TThermo,
    heNe2TThermo,
    pureZoneMixture,
    blottnerEucken2TTransport,
    sensibleInternalEnergy2T,
    eConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    ne2TThermo,
    heNe2TThermo,
    pureZoneMixture,
    blottnerEucken2TTransport,
    sensibleInternalEnergy2T,
    hConstThermo,
    perfectGas,
    specie
);

makeThermos
(
    ne2TThermo,
    heNe2TThermo,
    pureZoneMixture,
    blottnerEucken2TTransport,
    sensibleInternalEnergy2T,
    janaf2TThermo,
    perfectGas,
    specie
);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
