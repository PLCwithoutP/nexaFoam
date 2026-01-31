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

Description
    Kurganov-type flux calculation schemes.
    Kurganov-Tadmor and Kurganov-Noelle-Petrova schemes are given here.

\*---------------------------------------------------------------------------*/

//////////////////////////////////////////////////////////////////////
// --- 1) Directed interpolation of primitive fields onto faces --- //

surfaceScalarField rho_pos(interpolate(rho, pos)); // Full mixture
surfaceScalarField rho_neg(interpolate(rho, neg)); 

surfaceVectorField rhoU_pos(interpolate(rhoU, pos, U.name()));
surfaceVectorField rhoU_neg(interpolate(rhoU, neg, U.name()));

volScalarField rPsi("rPsi", 1.0/psi);
surfaceScalarField rPsi_pos(interpolate(rPsi, pos, TTR.name()));
surfaceScalarField rPsi_neg(interpolate(rPsi, neg, TTR.name()));

surfaceScalarField e_pos(interpolate(hTR, pos, TTR.name()));
surfaceScalarField e_neg(interpolate(hTR, neg, TTR.name()));

surfaceVectorField U_pos("U_pos", rhoU_pos/rho_pos);
surfaceVectorField U_neg("U_neg", rhoU_neg/rho_neg);

surfaceScalarField p_pos("p_pos", rho_pos*rPsi_pos);
surfaceScalarField p_neg("p_neg", rho_neg*rPsi_neg);

surfaceScalarField phi_fp("phi_fp", U_pos & mesh.Sf());
// Note: extracted out the orientation so becomes unoriented
phi_fp.setOriented(false);
surfaceScalarField phi_fn("phi_fn", U_neg & mesh.Sf());
phi_fn.setOriented(false);

// Make fluxes relative to mesh-motion
if (mesh.moving())
{
    surfaceScalarField meshPhi(mesh.phi());
    meshPhi.setOriented(false);
    phi_fp -= meshPhi;
    phi_fn -= meshPhi;
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// --- 2) Secondary calculations and interpolations --- //
volScalarField c("c", sqrt(thermo2T.CpTR()/thermo2T.CvT()*rPsi));
surfaceScalarField cSf_pos
(
    "cSf_pos",
    interpolate(c, pos, TTR.name())*mesh.magSf()
);

surfaceScalarField cSf_neg
(
    "cSf_neg",
    interpolate(c, neg, TTR.name())*mesh.magSf()
);

// Local speed of propagations
surfaceScalarField a_fp
(
    "a_fp",
    max(max(phi_fp + cSf_pos, phi_fn + cSf_neg), v_zero)
);

surfaceScalarField a_fn
(
    "a_fn",
    min(min(phi_fp - cSf_pos, phi_fn - cSf_neg), v_zero)
);

// Weighting factor, default KNP
surfaceScalarField alpha("alpha", a_fp/(a_fp - a_fn)); 

surfaceScalarField max_a("max_a", max(mag(a_fn), mag(a_fp)));

// Diffusive flux, default KNP
surfaceScalarField nw_f("nw_f", a_fn*alpha);

// KT corrections
if (fluxScheme == "Tadmor")
{
    alpha = 0.5; // Weighting factor
    nw_f = -0.5*max_a; // Diffusive flux
}
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// --- 3) Final flux calculation --- //
surfaceScalarField minus_alpha("minus_alpha", 1.0 - alpha);

// Weighted flux calculations
phi_fp *= alpha; // phi_f+ * alpha
phi_fn *= minus_alpha; // phi_f- * (1 - alpha)

// Diffusive corrections
surfaceScalarField aphi_fp("aphi_fp", phi_fp - nw_f);
surfaceScalarField aphi_fn("aphi_fn", phi_fn + nw_f);

// Reuse max_a for the maximum positive and negative fluxes
// estimated by the central scheme
// This value will be employed in delta t(timestep) calculations
max_a = max(mag(aphi_fp), mag(aphi_fn));
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// --- 4) Face fluxes for conservative variables --- //

// Mass Flux
phi = aphi_fp*rho_pos + aphi_fn*rho_neg;

// Convective Momentum Flux
surfaceVectorField phiU(aphi_fp*rhoU_pos + aphi_fn*rhoU_neg);
phiU.setOriented(true); // Note: reassembled orientation from the pos and neg parts so becomes oriented

// Full Momentum Flux
surfaceVectorField phiUp(phiU + (alpha*p_pos + minus_alpha*p_neg)*mesh.Sf());

// Total Energy Flux
surfaceScalarField phiEp
(
    "phiEp",
    aphi_fp*(rho_pos*(e_pos + 0.5*magSqr(U_pos)) + p_pos)
    + aphi_fn*(rho_neg*(e_neg + 0.5*magSqr(U_neg)) + p_neg)
    + nw_f*p_pos - nw_f*p_neg // pressure jump stabilization term
);

// Make flux for pressure-work absolute
if (mesh.moving())
{
    surfaceScalarField meshPhi(mesh.phi());
    meshPhi.setOriented(false);
    phiEp += meshPhi*(alpha*p_pos + minus_alpha*p_neg);
}

// Viscous Work Flux
surfaceScalarField sigmaDotU
(
    "sigmaDotU",
    (
        fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
        + fvc::dotInterpolate(mesh.Sf(), tauMC)
    )
    & (alpha*U_pos + minus_alpha*U_neg)
);
//////////////////////////////////////////////////////////////////////
