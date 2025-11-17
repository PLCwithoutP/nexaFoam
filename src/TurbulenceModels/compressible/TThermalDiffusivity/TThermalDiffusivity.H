/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2015-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::TThermalDiffusivity

Description
    Templated wrapper class to provide compressible turbulence models
    thermal diffusivity based thermal transport.

SourceFiles
    TThermalDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef TThermalDiffusivity_H
#define TThermalDiffusivity_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
               Class TThermalDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class TThermalDiffusivity
:
    public BasicTurbulenceModel
{

public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef volScalarField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    // Constructors

        //- Construct
        TThermalDiffusivity
        (
            const word& type,
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName
        );


    // Selectors

        //- Return a reference to the selected turbulence model
        static autoPtr<TThermalDiffusivity> New
        (
            const alphaField& alpha,
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


        //- Return a reference to the selected turbulence model
        static autoPtr<TThermalDiffusivity> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const surfaceScalarField& phi,
            const transportModel& transportModel,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~TThermalDiffusivity() = default;


    // Member Functions

        //- Return the thermal diffusivity for temperature [J/m/s/K]
        virtual tmp<volScalarField> kappa() const
        {
            return this->transport_.kappaTR();
        }

        //- Return the laminar thermal diffusivity for temperature on patch
        //  [J/m/s/K]
        virtual tmp<scalarField> kappa(const label patchi) const
        {
            return this->transport_.kappaTR(patchi);
        }

        //- Return the laminar thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alpha() const
        {
            return this->transport_.alpha();
        }

        //- Return the laminar thermal diffusivity for enthalpy on patch
        //  [kg/m/s]
        virtual tmp<scalarField> alpha(const label patchi) const
        {
            return this->transport_.alpha(patchi);
        }

        //- Return the turbulent thermal diffusivity for enthalpy [kg/m/s]
        virtual tmp<volScalarField> alphat() const;

        //- Return the turbulent thermal diffusivity for enthalpy for a patch
        //  [kg/m/s]
        virtual tmp<scalarField> alphat(const label patchi) const;

        //- Return the effective turbulent thermal diffusivity for temperature
        //  [J/m/s/K]
        virtual tmp<volScalarField> kappaEff() const
        {
            return kappa();
        }

        //- Return the effective turbulent thermal diffusivity for temperature
        //  [J/m/s/K]
        virtual tmp<scalarField> kappaEff(const label patchi) const
        {
            return kappa(patchi);
        }

        //- Return the effective turbulent thermal diffusivity for enthalpy
        //  [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
            return this->transport_.alphaheTR();
        }

        //- Return the effective turbulent thermal diffusivity for enthalpy
        //  for a patch [kg/m/s]
        virtual tmp<scalarField> alphaEff(const label patchi) const
        {
            return this->transport_.alphaheTR(patchi);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TThermalDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
