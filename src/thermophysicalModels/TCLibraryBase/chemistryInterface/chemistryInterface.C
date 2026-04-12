/*---------------------------------------------------------------------------*\
    chemistryInterface.C

    Factory for chemistryInterface::New().

    The dynamic_cast from tcLibraryInterface to ne2TReactionThermo lives
    here and nowhere else. When Mutation++ is implemented, an additional
    branch casts to the Mutation++ wrapper type instead.
\*---------------------------------------------------------------------------*/

#include "chemistryInterface.H"
#include "neTCChemistryAdapter.H"
#include "tcLibraryInterface.H"
#include "ne2TReactionThermo.H"

namespace Foam
{

autoPtr<chemistryInterface> chemistryInterface::New(tcLibraryInterface& thermo)
{
    ne2TReactionThermo& nativeThermo =
        dynamic_cast<ne2TReactionThermo&>(thermo);

    return autoPtr<chemistryInterface>
    (
        new neTCChemistryAdapter(nativeThermo)
    );
}

} // End namespace Foam