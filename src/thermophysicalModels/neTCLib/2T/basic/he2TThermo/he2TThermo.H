/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2021 OpenCFD Ltd.
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
    Foam::he2TThermo

Description
    Enthalpy/Internal energy for a mixture

SourceFiles
    he2TThermo.C

\*---------------------------------------------------------------------------*/

#ifndef he2TThermo_H
#define he2TThermo_H

#include "basicMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class he2TThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Basic2TThermo, class MixtureType>
class he2TThermo
:
    public Basic2TThermo,
    public MixtureType
{
protected:

    // Protected data

        //- Energy field
        volScalarField h_;


    // Protected Member Functions

        // Enthalpy/Internal energy

            //- Correct the enthalpy/internal energy field boundaries
            void hBoundaryCorrection(volScalarField& h);


private:

    // Private Member Functions

        //- Construct as copy (not implemented)
        he2TThermo(const he2TThermo<Basic2TThermo, MixtureType>&);


        //- Initialize he2TThermo
        void init
        (
            const volScalarField& p,
            const volScalarField& TTR,
            const volScalarField& TVib,
            volScalarField& h
        );


public:

    // Constructors

        //- Construct from mesh
        he2TThermo
        (
            const fvMesh&,
            const word& phaseName
        );

        //- Construct from mesh and dictionary
        he2TThermo
        (
            const fvMesh&,
            const dictionary&,
            const word& phaseName
        );

        //- Construct from mesh,dictionary,phase name with a single temperature
        he2TThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictionaryName
        );


    //- Destructor
    virtual ~he2TThermo();


    // Member functions

        //- Return the composition of the mixture
        virtual typename MixtureType::basicMixtureType&
        composition()
        {
            return *this;
        }

        //- Return the composition of the mixture
        virtual const typename MixtureType::basicMixtureType&
        composition() const
        {
            return *this;
        }

        //- Return the name of the thermo physics
        virtual word thermoName() const
        {
            return MixtureType::thermoType::typeName();
        }

        //- Return true if the equation of state is incompressible
        //  i.e. rho != f(p)
        virtual bool incompressible() const
        {
            return MixtureType::thermoType::incompressible;
        }

        //- Return true if the equation of state is isochoric
        //  i.e. rho = const
        virtual bool isochoric() const
        {
            return MixtureType::thermoType::isochoric;
        }


        // Access to thermodynamic state variables

            //- Enthalpy/Internal energy [J/kg]
            //  Non-const access allowed for transport equations
            virtual volScalarField& h()
            {
                return h_;
            }

            //- Enthalpy/Internal energy [J/kg]
            virtual const volScalarField& h() const
            {
                return h_;
            }


        // Fields derived from thermodynamic state variables

            //- Enthalpy/Internal energy
            //  for given pressure and temperature [J/kg]
            virtual tmp<volScalarField> h
            (
                const volScalarField& p,
                const volScalarField& TTR
            ) const;

            //- Enthalpy/Internal energy for cell-set [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Enthalpy/Internal energy for patch [J/kg]
            virtual tmp<scalarField> h
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Chemical enthalpy [J/kg]
            virtual tmp<volScalarField> hc() const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> THE_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> THE_TR
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const label patchi
            ) const;

            //- Temperature from enthalpy/internal energy for cell-set
            virtual tmp<scalarField> TEH_Vib
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const labelList& cells
            ) const;

            //- Temperature from enthalpy/internal energy for patch
            virtual tmp<scalarField> TEH_Vib
            (
                const scalarField& h,
                const scalarField& p,
                const scalarField& TTR0,      // starting translational-rotational temperature
                const scalarField& TVib0,      // starting vibrational temperature
                const label patchi
            ) const;

            //- Heat capacity at constant pressure for patch [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Heat capacity using pressure and temperature [J/kg/K]
            virtual tmp<scalarField> CpTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Heat capacity at constant pressure [J/kg/K]
            virtual tmp<volScalarField> CpTR() const;

            //- Translational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvT() const;

            //- Translational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvT
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Rotational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvR() const;

            //- Rotational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Vibrational heat capacity at constant volume [J/kg/K]
            virtual tmp<volScalarField> CvVib() const;

            //- Vibrational heat capacity at constant volume for patch [J/kg/K]
            virtual tmp<scalarField> CvVib
            (
                const scalarField& p,
                const scalarField& TTR,
                const scalarField& Tvib,
                const scalar ThetaVib,
                const label patchi
            ) const;

            //- Density from pressure and temperature
            virtual tmp<scalarField> rhoEoS
            (
                const scalarField& p,
                const scalarField& TTR,
                const labelList& cells
            ) const;

            //- Gamma = CpTR/CvT []
            virtual tmp<volScalarField> gamma() const;

            //- Gamma = CpTR/CvT for patch []
            virtual tmp<scalarField> gamma
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Translational-Rotational heat capacity at constant pressure/volume for patch [J/kg/K]
            virtual tmp<scalarField> CpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Translational-Rotational heat capacity at constant pressure/volume [J/kg/K]
            virtual tmp<volScalarField> CpvTR() const;

            //- Translational-Rotational heat capacity ratio []
            virtual tmp<volScalarField> CpByCpvTR() const;

            //- Translational-Rotational heat capacity ratio for patch []
            virtual tmp<scalarField> CpByCpvTR
            (
                const scalarField& p,
                const scalarField& TTR,
                const label patchi
            ) const;

            //- Molecular weight [kg/kmol]
            virtual tmp<volScalarField> W() const;


        // Fields derived from transport state variables

            //- Translational-rotational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaTR() const;

            //- Translational-rotational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaTR
            (
                const label patchi
            ) const;

            //- Vibrational thermal diffusivity for temperature of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaVib() const;

            //- Vibrational thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaVib
            (
                const label patchi
            ) const;

             //- Translational-rotational thermal diffusivity for energy of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaheTR() const;

            //- Translational-rotational thermal diffusivity for energy of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaheTR(const label patchi) const;

            //- Effective thermal diffusivity for temperature
            //  of mixture [J/m/s/K]
            virtual tmp<volScalarField> kappaEff(const volScalarField&) const;

            //- Effective thermal diffusivity for temperature
            //  of mixture for patch [J/m/s/K]
            virtual tmp<scalarField> kappaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;

            //- Effective thermal diffusivity of mixture [kg/m/s]
            virtual tmp<volScalarField> alphaEff
            (
                const volScalarField& alphat
            ) const;

            //- Effective thermal diffusivity of mixture for patch [kg/m/s]
            virtual tmp<scalarField> alphaEff
            (
                const scalarField& alphat,
                const label patchi
            ) const;


        //- Read thermophysical properties dictionary
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

#ifdef NoRepository
    #include "he2TThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
