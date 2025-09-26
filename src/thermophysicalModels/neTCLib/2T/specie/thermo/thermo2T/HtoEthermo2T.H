inline scalar CvT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar cpTR = CpTR(p, TTR);
    volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
    return cpTR - cpmcvTR;
    #else
    return CpTR(p, TTR) - EquationOfState::CpMCv(p, TTR);
    #endif
}

inline scalar CvR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR);
    }
    #endif
}

inline scalar CvVib
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar dimensionlessVibTemp = 0.5*ThetaVib/TVib 
        volatile const scalar div = sinh(0.5*ThetaVib/TVib);
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR*dimensionlessVibTemp/div;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR)*(0.5*ThetaVib/TVib)/sinh(0.5*ThetaVib/TVib);
    }
    #endif
}

inline scalar EsT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar hsT = Hs(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return hsT - p/rho;
    #else
    return Hs(p, TTR) - p/EquationOfState::rho(p, TTR);
    #endif
}

inline scalar EsR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR*TTR;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR)*TTR;
    }
    #endif
}

inline scalar EsVib
(
    const scalar p,
    const scalar TTR,
    const scalar TVib,
    const scalar ThetaVib
) const
{
    #ifdef __clang__
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        volatile const scalar div = exp(ThetaVib/TVib) - 1;
        volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
        return cpmcvTR*ThetaVib/div;
    }
    #else
    if (static_cast<const specie2T&>(*this).isAtomic())
    {
        return 0;
    }
    else
    {
        return EquationOfState::CpMCv(p, TTR)*ThetaVib/(exp(ThetaVib/TVib) - 1);
    }
    #endif
}

inline scalar EaT
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar haT = Ha(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return haT - p/rho;
    #else
    return Ha(p, TTR) - p/EquationOfState::rho(p, TTR);
    #endif
}
