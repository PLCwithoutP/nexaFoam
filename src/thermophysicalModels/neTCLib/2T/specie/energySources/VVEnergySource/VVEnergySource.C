/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "VVEnergySource.H"
#include "fvMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class MixtureType, class MixingRule>
Foam::VVEnergySource<MixtureType, MixingRule>::VVEnergySource
(
    MixtureType& mixture,
    MixingRule& MR
)
:
    baseEnergySource<MixtureType, MixingRule>(mixture,MR),
    mix_(mixture),
    mr_(MR),
    names_(mixture.species()),
    sigma_prime_ (3e-21), // N2, O2, NO
    boltzmann_const_ (1.380649e-23),
    avagadro_const_ (6.022e23),
    Pml_const_ (1e-2),
    Q_VVList_()
{
     
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class MixtureType, class MixingRuleType>
void Foam::VVEnergySource<MixtureType, MixingRuleType>::makeQVibSourceFields
(
    const fvMesh& mesh
)
{
    if (Q_VVList_.size())
    {
        // already built
        return;
    }

    const PtrList<volScalarField>& Y_species = mix_.Y();
    Q_VVList_.setSize(Y_species.size());

    forAll(Y_species, i)
    {
        const volScalarField& Yi = Y_species[i];

        const word fieldName("Q_VV_" + Yi.name());

        Q_VVList_.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    fieldName,
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar
                (
                    "zero",
                    dimEnergy/dimVolume/dimTime,
                    0.0
                )
            )
        );
    }
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::sigmai
(
    const scalar TTR 
)
{
    //Info << "Cross Section: " << sigma_prime_*(50000.0/TTR)*(50000.0/TTR) << nl; 
    return sigma_prime_*(50000.0/TTR)*(50000.0/TTR);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::c_bar_s
(
    const scalar TTR,
    const label s 
)
{
    //Info << "Most Probable Speed: " << pow((8*mix_.R(s)*TTR)/(constant::mathematical::pi), 0.5) << nl; 
    return pow((8*mix_.R(s)*TTR)/(constant::mathematical::pi) , 0.5);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::n_total
(
    const scalar p,
    const scalar TTR
)
{
    //Info << "Total Number Density: " << p/(TTR*boltzmann_const_) << nl; 
    return p/(TTR*boltzmann_const_);
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::n_s
(
    const scalar p,
    const scalar TTR,
    const label s,
    const label celli
)
{
    // total number density
    const scalar nTot = n_total(p, TTR);   // calls the 2-argument version

    // species mole fraction in this cell
    volScalarField& XrField = mr_.computeXiFromYi(s);
    const scalar Xr = XrField[celli];

    // species number density
    const scalar nSpecie = Xr*nTot;

    /* Info<< "Number density of specie " << s
        << " at cell " << celli
        << " = " << nSpecie
        << " [1/m3], Xr = " << Xr
        << ", n_tot = " << nTot << nl; */

    return nSpecie;
}

template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::p_s
(
    const scalar p,
    const scalar TTR,   // kept for symmetry; not used
    const label s,
    const label celli
)
{
    volScalarField& XsField = mr_.computeXiFromYi(s);
    const scalar Xs = XsField[celli];

    // partial pressure: p_s = X_s * p_total
    const scalar pSpecie = Xs*p;

    /* Optional debug:
    Info<< "Partial pressure of specie " << s
        << " at cell " << celli
        << " = " << pSpecie
        << " Pa (Xs=" << Xs << ", p_tot=" << p << ')' << nl;
    */

    return pSpecie;
}


// Knab's V–V formula: Q_{m,V-V} for one species m (=s)
template<class MixtureType, class MixingRule>
Foam::scalar 
Foam::VVEnergySource<MixtureType, MixingRule>::Q_VV_s
(
    const scalar p,      
    const scalar TTR,
    const scalar TVib,
    const label  s,      
    const label  celli
)
{
    scalar Qm = 0.0;

    const label nSpec = mix_.Y().size();

    const scalar p_m   = p_s(p, TTR, s, celli);          
    const scalar rho_m = mix_.rho(s, p_m, TTR);          

    const scalar ev_s_Ttr =
        mix_.EsVib(s, p_m, TTR, TTR, thetai(s));         
    const scalar ev_s_Tvm =
        mix_.EsVib(s, p_m, TTR, TVib, thetai(s));        

    const scalar Mm = Wi(s)*1e-3;                             
    const scalar cbar_m = c_bar_s(TTR, s);               
    Info << "Most probable speed is : " << cbar_m << nl;
    
    for (label r = 0; r < nSpec; ++r)
    {
        if (r == s) continue;                            // r ≠ s
        if (!mix_.isSpecieMolecular(r)) continue;

        const scalar Ml = Wi(r)*1e-3;                         

        const scalar p_r   = p_s(p, TTR, r, celli);      
        const scalar rho_l = mix_.rho(r, p_r, TTR);      


        const scalar velFactor =
            cbar_m * Foam::sqrt(rho_l/Ml);

        Info << "Vel factor is : " << velFactor << nl;

        //const scalar sigma_ml = sigmai(TTR);
        const scalar sigma_ml = sigma_prime_;            
        const scalar P_ml     = Pml_const_;              

        const scalar ev_r_Ttr =
            mix_.EsVib(r, p_r, TTR, TTR, thetai(r));     
        const scalar ev_r_Tvl =
            mix_.EsVib(r, p_r, TTR, TVib, thetai(r));    
        const scalar energyBracket =
            ev_s_Ttr*(ev_r_Tvl/ev_r_Ttr) - ev_s_Tvm;

        Info << "Energy bracket is : " << energyBracket << nl;

        Qm +=
            avagadro_const_   
          * sigma_ml          
          * P_ml              
          * velFactor
          * rho_m
          * energyBracket;
    }

    return Qm;   
}


template<class MixtureType, class MixingRule>
Foam::PtrList<Foam::volScalarField>&
Foam::VVEnergySource<MixtureType, MixingRule>::correctVibVibSource
(
    const volScalarField& p,
    const volScalarField& TTR,
    const volScalarField& TVib
)
{
    makeQVibSourceFields(TTR.mesh());

    const scalarField& pCells    = p.primitiveField();
    const scalarField& TTRCells  = TTR.primitiveField();
    const scalarField& TVibCells = TVib.primitiveField();

    const PtrList<volScalarField>& Y_species = mix_.Y();

    forAll(Y_species, m)
    {
        scalarField& Q_VV_Cells = Q_VVList_[m].primitiveFieldRef();

        forAll(TTRCells, celli)
        {
            const scalar pCell    = pCells[celli];
            const scalar TTRCell  = TTRCells[celli];
            const scalar TVibCell = TVibCells[celli];

            Q_VV_Cells[celli] = Q_VV_s
            (
                pCell,      
                TTRCell,
                TVibCell,
                m,
                celli
            );
        }

        Q_VVList_[m].correctBoundaryConditions();
    }

    return Q_VVList_;
}


// ************************************************************************* //
