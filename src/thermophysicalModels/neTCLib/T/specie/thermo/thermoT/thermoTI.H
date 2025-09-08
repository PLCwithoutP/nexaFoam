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

#include "thermoT.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
inline Foam::species::thermoT<ThermoT, Type>::thermoT
(
    const ThermoT& sp
)
:
    ThermoT(sp)
{}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (thermoT<ThermoT, Type>::*F)(const scalar, const scalar) const,
    scalar (thermoT<ThermoT, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (thermoT<ThermoT, Type>::*limit)(const scalar) const
) const
{
    if (T0 < 0)
    {
        FatalErrorInFunction
            << "Negative initial temperature T0: " << T0
            << abort(FatalError);
    }

    scalar Test = T0;
    scalar Tnew = T0;
    scalar Ttol = T0*tol_;
    int    iter = 0;

    do
    {
        Test = Tnew;
        Tnew =
            (this->*limit)
            (Test - ((this->*F)(p, Test) - f)/(this->*dFdT)(p, Test));

        if (iter++ > maxIter_)
        {
            FatalErrorInFunction
                << "Maximum number of iterations exceeded: " << maxIter_
                << " when starting from T0:" << T0
                << " old T:" << Test << " new T:" << Tnew
                << " f:" << f
                << " p:" << p
                << " tol:" << Ttol
                << abort(FatalError);
        }

    } while (mag(Tnew - Test) > Ttol);

    return Tnew;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
inline Foam::species::thermoT<ThermoT, Type>::thermoT
(
    const word& name,
    const thermoT& st
)
:
    ThermoT(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
inline Foam::word
Foam::species::thermoT<ThermoT, Type>::heName()
{
    return Type<thermoT<ThermoT, Type>>::energyName();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<thermoT<ThermoT, Type>>::Cpv(*this, p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::gamma(const scalar p, const scalar T) const
{
    #ifdef __clang__
    volatile const scalar Cp = this->Cp(p, T);
    #else
    const scalar Cp = this->Cp(p, T);
    #endif

    return Cp/(Cp - this->CpMCv(p, T));
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::CpByCpv
(
    const scalar p,
    const scalar T
) const
{
    return Type<thermoT<ThermoT, Type>>::CpByCpv(*this, p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::HE(const scalar p, const scalar T) const
{
    return Type<thermoT<ThermoT, Type>>::HE(*this, p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T*this->S(p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T*this->S(p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T)*this->W();
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::K(const scalar p, const scalar T) const
{
    scalar arg = -this->Y()*this->G(Pstd, T)/(RR*T);

    if (arg < 600)
    {
        return exp(arg);
    }
    else
    {
        return VGREAT;
    }
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::Kc(const scalar p, const scalar T) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/(RR*T), nm);
    }
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::Kx
(
    const scalar p,
    const scalar T
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(Pstd/p, nm);
    }
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::Kn
(
    const scalar p,
    const scalar T,
    const scalar n
) const
{
    const scalar nm = this->Y()/this->W();

    if (equal(nm, SMALL))
    {
        return Kp(p, T);
    }
    else
    {
        return Kp(p, T)*pow(n*Pstd/p, nm);
    }
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    return Type<thermoT<ThermoT, Type>>::THE(*this, he, p, T0);
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::THs
(
    const scalar hs,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        hs,
        p,
        T0,
        &thermoT<ThermoT, Type>::Hs,
        &thermoT<ThermoT, Type>::Cp,
        &thermoT<ThermoT, Type>::limit
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::THa
(
    const scalar ha,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ha,
        p,
        T0,
        &thermoT<ThermoT, Type>::Ha,
        &thermoT<ThermoT, Type>::Cp,
        &thermoT<ThermoT, Type>::limit
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::TEs
(
    const scalar es,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        es,
        p,
        T0,
        &thermoT<ThermoT, Type>::Es,
        &thermoT<ThermoT, Type>::Cv,
        &thermoT<ThermoT, Type>::limit
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar Foam::species::thermoT<ThermoT, Type>::TEa
(thermoT_H
    const scalar ea,
    const scalar p,
    const scalar T0
) const
{
    return T
    (
        ea,
        p,
        T0,
        &thermoT<ThermoT, Type>::Ea,
        &thermoT<ThermoT, Type>::Cv,
        &thermoT<ThermoT, Type>::limit
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::dKcdTbyKc
(
    const scalar p,
    const scalar T
) const
{
    const scalar dKcdTbyKc =
        (this->S(Pstd, T) + this->Gstd(T)/T)*this->Y()/(RR*T);

    const scalar nm = this->Y()/this->W();
    if (equal(nm, SMALL))
    {
        return dKcdTbyKc;
    }
    else
    {
        return dKcdTbyKc - nm/T;
    }
}


template<class ThermoT, template<class> class Type>
inline Foam::scalar
Foam::species::thermoT<ThermoT, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
inline void Foam::species::thermoT<ThermoT, Type>::operator+=
(
    const thermoT<ThermoT, Type>& st
)
{
    ThermoT::operator+=(st);
}


template<class ThermoT, template<class> class Type>
inline void Foam::species::thermoT<ThermoT, Type>::operator*=(const scalar s)
{
    ThermoT::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class ThermoT, template<class> class Type>
inline Foam::species::thermoT<ThermoT, Type> Foam::species::operator+
(
    const thermoT<ThermoT, Type>& st1,
    const thermoT<ThermoT, Type>& st2
)
{
    return thermoT<ThermoT, Type>
    (
        static_cast<const ThermoT&>(st1) + static_cast<const ThermoT&>(st2)
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::species::thermoT<ThermoT, Type> Foam::species::operator*
(
    const scalar s,
    const thermoT<ThermoT, Type>& st
)
{
    return thermoT<ThermoT, Type>
    (
        s*static_cast<const ThermoT&>(st)
    );
}


template<class ThermoT, template<class> class Type>
inline Foam::species::thermoT<ThermoT, Type> Foam::species::operator==
(
    const thermoT<ThermoT, Type>& st1,
    const thermoT<ThermoT, Type>& st2
)
{
    return thermoT<ThermoT, Type>
    (
        static_cast<const ThermoT&>(st1) == static_cast<const ThermoT&>(st2)
    );
}


// ************************************************************************* //
