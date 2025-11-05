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

\*---------------------------------------------------------------------------*/

#include "WilkeMR.H"
#include "volFields.H"
#include "fvc.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

// --- local utilities --------------------------------------------------------

// Put near the top of WilkeMR.C (local to this TU)
static inline scalar wilke_phi_scalar
(
    const scalar mu_i,  // species i dynamic viscosity
    const scalar mu_j,  // species j dynamic viscosity
    const scalar Mi,    // species i molecular weight
    const scalar Mj     // species j molecular weight
)
{
    // Guard tiny/zero values for stability
    const scalar mujSafe = max(SMALL, mu_j);

    // Classic Wilke φ_ij:
    // φ_ij = [1 + sqrt(μ_i/μ_j) * (M_j/M_i)^(1/4)]^2
    //        / [ sqrt(8) * sqrt(1 + M_i/M_j) ]
    const scalar muRatioSqrt = std::sqrt(mu_i/mujSafe);
    const scalar Mi_over_Mj  = Mi/Mj;
    const scalar Mj_over_Mi  = Mj/Mi;

    const scalar term = 1.0 + muRatioSqrt * std::pow(Mj_over_Mi, 0.25);
    return sqr(term) / (std::sqrt(8.0) * std::sqrt(1.0 + Mi_over_Mj));
}

// --- WilkeMR methods --------------------------------------------------------

WilkeMR::WilkeMR
(
    const basicMultiComponent2TMixture& mc2Tmix
)
:
    baseMixingRule(mc2Tmix),
    W_(),
    mu_(),
    x_()
{
    // Pull species mass fractions to infer species count & names
    const PtrList<volScalarField>& Y = mixture_.Y();  // expected API
    const label nSpec = Y.size();

    // sizes
    W_.setSize(nSpec);
    mu_.setSize(nSpec);
    x_.setSize(nSpec);

    // molecular weights (assumed API: mixture_.W(i) -> scalar)
    for (label i=0; i<nSpec; ++i)
    {
        W_[i] = mixture_.W(i);
    }

    // species viscosities: assumed API:
    //   - either mixture_.mu(i) -> volScalarField,
    //   - or mixture_.speciesMu(i) -> volScalarField.
    // If your API differs, please tell me the correct accessor.
    for (label i=0; i<nSpec; ++i)
    {
        mu_.set
        (
            i,
            new volScalarField
            (
                // try mixture_.mu(i); change here if your accessor differs
                mixture_.mu(i)
            )
        );
    }

    // initialize x_ from current Y
    updateX();
}

WilkeMR::~WilkeMR() 
{

}

// find index of Yi by pointer identity
label WilkeMR::specieIndex
(
    const volScalarField& Yi
) const
{
    const PtrList<volScalarField>& Y = mixture_.Y();
    for (label i=0; i<Y.size(); ++i)
    {
        if (&Y[i] == &Yi) return i;
    }
    FatalErrorInFunction
        << "Given Yi is not one of mixture_.Y() fields." << nl
        << abort(FatalError);
    return -1;
}

// recompute x_k from current Y_k
void WilkeMR::updateX
()
{
    const PtrList<volScalarField>& Y = mixture_.Y();
    const label nSpec = Y.size();

    // denom = sum_k (Y_k / W_k)
    tmp<volScalarField> tDen
    (
        volScalarField::New
        (
            "sumYoverW",
            Y[0].mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    volScalarField& den = tDen.ref();

    for (label k=0; k<nSpec; ++k)
    {
        den += Y[k] / W_[k];
    }

    for (label k=0; k<nSpec; ++k)
    {
        x_.set
        (
            k,
            new volScalarField
            (
                (Y[k]/W_[k]) / den
            )
        );
        x_[k].rename("x_"+Y[k].name());
    }
}

// φ_ij (field)
tmp<volScalarField> WilkeMR::phi_ij
(
    const label i, 
    const label j
) const
{
    const volScalarField& mui = mu_[i];
    const volScalarField& muj = mu_[j];

    const scalar Mi = W_[i];
    const scalar Mj = W_[j];

    // sqrt(μ_i/μ_j) – field quantity
    tmp<volScalarField> tRatio = sqrt(mui/muj);

    // scalar factors with mol. weights
    const scalar Mi_over_Mj = Mi/Mj;
    const scalar Mj_over_Mi = Mj/Mi;

    // Build φ_ij field
    tmp<volScalarField> tPhi
    (
        volScalarField::New
        (
            "phi_"+name(i)+"_"+name(j),
            mui.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    volScalarField& phi = tPhi.ref();

    // compute cellwise using scalar weight factors and field sqrt(mu_i/mu_j)
    phi = sqr(1.0 + tRatio() * pow(Mj_over_Mi, 0.25))
        / (sqrt(8.0)*sqrt(1.0 + Mi_over_Mj));

    return tPhi;
}

// Y -> X for the *same* species
tmp<volScalarField> WilkeMR::YtoX
(
    const volScalarField& Yi
) const
{
    const label i = specieIndex(Yi);

    const PtrList<volScalarField>& Y = mixture_.Y();
    const label nSpec = Y.size();

    tmp<volScalarField> tDen
    (
        volScalarField::New
        (
            "sumYoverW",
            Yi.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    volScalarField& den = tDen.ref();

    for (label k=0; k<nSpec; ++k)
    {
        den += Y[k] / W_[k];
    }

    tmp<volScalarField> tXi
    (
        volScalarField::New
        (
            "x_"+Yi.name(),
            Yi.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );

    tXi.ref() = (Yi/W_[i]) / den;
    return tXi;
}

tmp<volScalarField> WilkeMR::Qscalar
(
    const volScalarField& Qspecie,   // Q_s field (e.g. mu_i)
    const volScalarField& Yi         // identifies specie s
) const
{
    const label i = specieIndex(Yi);

    // keep mole fractions current
    const_cast<WilkeMR*>(this)->updateX();

    const PtrList<volScalarField>& Y = mixture_.Y();
    const label nSpec = Y.size();

    // φ_i = X_i + Σ_{j≠i} X_j Φ_ij   (dimensionless)
    tmp<volScalarField> tDen
    (
        volScalarField::New
        (
            "wilkeDen_"+name(i),
            Yi.mesh(),
            dimensionedScalar("zero", dimless, 0.0)
        )
    );
    volScalarField& den = tDen.ref();
    den = x_[i];  // X_i term

    const volScalarField& mui = mu_[i];
    const scalar Mi = W_[i];

    for (label j=0; j<nSpec; ++j)
    {
        if (j == i) continue;

        const volScalarField& muj = mu_[j];
        const scalar Mj = W_[j];

        const dimensionedScalar tinyMu("tinyMu", muj.dimensions(), SMALL);

        tmp<volScalarField> tPhi
        (
            volScalarField::New
            (
                "Phi_"+name(i)+"_"+name(j),
                Yi.mesh(),
                dimensionedScalar("zero", dimless, 0.0)
            )
        );

        // Φ_ij (field):
        tPhi.ref() =
            sqr( 1.0 + sqrt(mui / max(muj, tinyMu)) * pow(Mj/Mi, 0.25) )
          / ( sqrt(8.0) * sqrt(1.0 + Mi/Mj) );

        den += x_[j] * tPhi;
    }

    // Return Q_s * X_s / φ_s (dimensions = Q)
    tmp<volScalarField> tOut
    (
        volScalarField::New
        (
            "Qmix_contrib_"+Yi.name(),
            Yi.mesh(),
            Qspecie.dimensions()  // ensure result has Q’s dimensions
        )
    );

    tOut.ref() = ( x_[i] * Qspecie )
               / max(den, dimensionedScalar("epsDen", dimless, SMALL));

    return tOut;
}



tmp<scalarField> WilkeMR::Qscalar
(
    const label patchi,
    const volScalarField& Qspecie,   // Q_s field
    const volScalarField& Yi
) const
{
    const label i = specieIndex(Yi);

    const_cast<WilkeMR*>(this)->updateX();

    const PtrList<volScalarField>& Y = mixture_.Y();
    const label nSpec = Y.size();

    const scalarField& Xi  = x_[i].boundaryField()[patchi];
    const scalarField& mui = mu_[i].boundaryField()[patchi];
    const scalarField& Qi  = Qspecie.boundaryField()[patchi];
    const scalar Mi = W_[i];

    // φ_i(face) = X_i + Σ_{j≠i} X_j Φ_ij
    scalarField den = Xi;

    for (label j=0; j<nSpec; ++j)
    {
        if (j == i) continue;

        const scalarField& xj  = x_[j].boundaryField()[patchi];
        const scalarField& muj = mu_[j].boundaryField()[patchi];
        const scalar Mj = W_[j];

        scalarField Phi(mui.size());
        forAll(Phi, c)
        {
            Phi[c] = wilke_phi_scalar(
                mui[c],  // μ_i at face c
                muj[c],  // μ_j at face c
                Mi,      // M_i
                Mj       // M_j
            );
        }

        den += xj * Phi;
    }

    // Return Q_s * X_s / φ_s per face
    tmp<scalarField> tOut(new scalarField(Qi.size()));
    scalarField& out = tOut.ref();

    forAll(out, c)
    {
        const scalar d = max(den[c], VSMALL);
        out[c] = Xi[c]*Qi[c] / d;
    }

    return tOut;
}


} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::WilkeMR::WilkeMR
()
:
    W_ (),
    mu_ (),
    x_ ()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::WilkeMR::~WilkeMR()
{}


// ************************************************************************* //