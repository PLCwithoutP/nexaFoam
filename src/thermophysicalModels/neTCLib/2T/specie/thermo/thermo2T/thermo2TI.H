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

#include "thermo2T.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type>::thermo2T
(
    const Thermo2T& sp
)
:
    Thermo2T(sp)
{}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::T
(
    scalar f,
    scalar p,
    scalar T0,
    scalar (thermo2T<Thermo2T, Type>::*F)(const scalar, const scalar) const,
    scalar (thermo2T<Thermo2T, Type>::*dFdT)(const scalar, const scalar)
        const,
    scalar (thermo2T<Thermo2T, Type>::*limit)(const scalar) const
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

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type>::thermo2T
(
    const word& name,
    const thermo2T& st
)
:
    Thermo2T(name, st)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::word
Foam::species::thermo2T<Thermo2T, Type>::heName()
{
    return Type<thermo2T<Thermo2T, Type>>::energyName();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::Cpv(const scalar p, const scalar T) const
{
    return Type<thermo2T<Thermo2T, Type>>::Cpv(*this, p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::gamma(const scalar p, const scalar T) const
{
    #ifdef __clang__
    volatile const scalar Cp = this->Cp(p, T);
    #else
    const scalar Cp = this->Cp(p, T);
    #endif

    return Cp/(Cp - this->CpMCv(p, T));
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::CpByCpv
(
    const scalar p,
    const scalar T
) const
{
    return Type<thermo2T<Thermo2T, Type>>::CpByCpv(*this, p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::HE(const scalar p, const scalar T) const
{
    return Type<thermo2T<Thermo2T, Type>>::HE(*this, p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::G(const scalar p, const scalar T) const
{
    return this->Ha(p, T) - T*this->S(p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::A(const scalar p, const scalar T) const
{
    return this->Ea(p, T) - T*this->S(p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cp(const scalar p, const scalar T) const
{
    return this->Cp(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::ha(const scalar p, const scalar T) const
{
    return this->Ha(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::hs(const scalar p, const scalar T) const
{
    return this->Hs(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::hc() const
{
    return this->Hc()*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::s(const scalar p, const scalar T) const
{
    return this->S(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::he(const scalar p, const scalar T) const
{
    return this->HE(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::cv(const scalar p, const scalar T) const
{
    return this->Cv(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::es(const scalar p, const scalar T) const
{
    return this->Es(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::ea(const scalar p, const scalar T) const
{
    return this->Ea(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::g(const scalar p, const scalar T) const
{
    return this->G(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::a(const scalar p, const scalar T) const
{
    return this->A(p, T)*this->W();
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::K(const scalar p, const scalar T) const
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


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::Kp(const scalar p, const scalar T) const
{
    return K(p, T);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::Kc(const scalar p, const scalar T) const
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


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::Kx
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


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::Kn
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


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::THE
(
    const scalar he,
    const scalar p,
    const scalar T0
) const
{
    return Type<thermo2T<Thermo2T, Type>>::THE(*this, he, p, T0);
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::THs
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
        &thermo2T<Thermo2T, Type>::Hs,
        &thermo2T<Thermo2T, Type>::Cp,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::THa
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
        &thermo2T<Thermo2T, Type>::Ha,
        &thermo2T<Thermo2T, Type>::Cp,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TEs
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
        &thermo2T<Thermo2T, Type>::Es,
        &thermo2T<Thermo2T, Type>::Cv,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar Foam::species::thermo2T<Thermo2T, Type>::TEa
(thermo2T_H
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
        &thermo2T<Thermo2T, Type>::Ea,
        &thermo2T<Thermo2T, Type>::Cv,
        &thermo2T<Thermo2T, Type>::limit
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::dKcdTbyKc
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


template<class Thermo2T, template<class> class Type>
inline Foam::scalar
Foam::species::thermo2T<Thermo2T, Type>::dcpdT(const scalar p, const scalar T) const
{
    return this->dCpdT(p, T)*this->W();
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline void Foam::species::thermo2T<Thermo2T, Type>::operator+=
(
    const thermo2T<Thermo2T, Type>& st
)
{
    Thermo2T::operator+=(st);
}


template<class Thermo2T, template<class> class Type>
inline void Foam::species::thermo2T<Thermo2T, Type>::operator*=(const scalar s)
{
    Thermo2T::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator+
(
    const thermo2T<Thermo2T, Type>& st1,
    const thermo2T<Thermo2T, Type>& st2
)
{
    return thermo2T<Thermo2T, Type>
    (
        static_cast<const Thermo2T&>(st1) + static_cast<const Thermo2T&>(st2)
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator*
(
    const scalar s,
    const thermo2T<Thermo2T, Type>& st
)
{
    return thermo2T<Thermo2T, Type>
    (
        s*static_cast<const Thermo2T&>(st)
    );
}


template<class Thermo2T, template<class> class Type>
inline Foam::species::thermo2T<Thermo2T, Type> Foam::species::operator==
(
    const thermo2T<Thermo2T, Type>& st1,
    const thermo2T<Thermo2T, Type>& st2
)
{
    return thermo2T<Thermo2T, Type>
    (
        static_cast<const Thermo2T&>(st1) == static_cast<const Thermo2T&>(st2)
    );
}


// ************************************************************************* //
