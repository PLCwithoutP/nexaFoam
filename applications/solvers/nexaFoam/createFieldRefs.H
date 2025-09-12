volScalarField& p = thermo2T.p();
const volScalarField& TTR = thermo2T.TTR();
const volScalarField& TVib = thermo2T.TVib();
const volScalarField& psi = thermo2T.psi();

bool inviscid(true);
if (max(thermo2T.mu().cref().primitiveField()) > 0.0)
{
    inviscid = false;
}
