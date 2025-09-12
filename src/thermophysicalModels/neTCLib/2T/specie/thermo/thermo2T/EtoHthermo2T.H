inline scalar CpTR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar cvT = CvT(p, TTR);
    volatile const scalar cpmcvTR = EquationOfState::CpMCv(p, TTR);
    return cvT + cpmcvTR;
    #else
    return CvT(p, TTR) + EquationOfState::CpMCv(p, TTR);
    #endif
}

inline scalar HsTR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar esTR = EsTR(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return esTR + p/rho;
    #else
    return EsTR(p, TTR) + p/EquationOfState::rho(p, TTR);
    #endif
}

inline scalar HaTR
(
    const scalar p,
    const scalar TTR
) const
{
    #ifdef __clang__
    volatile const scalar eaTR = EaTR(p, TTR);
    volatile const scalar rho = EquationOfState::rho(p, TTR);
    return eaTR + p/rho;
    #else
    return EaTR(p, TTR) + p/EquationOfState::rho(p, TTR);
    #endif
}