/*---------------------------------------------------------------------------*\
    tcLibraryInterface.C

    Factory implementation for tcLibraryInterface::New().

    Reads 'library' keyword from constant/thermophysicalProperties:
      library  neTCLib;     -> constructs ne2TReactionThermo
      library  default;     -> constructs ne2TReactionThermo (alias)
      library  mutationpp;  -> FatalError (not yet implemented)
      (keyword absent)      -> constructs ne2TReactionThermo (safe default)
\*---------------------------------------------------------------------------*/

#include "tcLibraryInterface.H"
#include "ne2TReactionThermo.H"
#include "IOdictionary.H"
#include "Time.H"

namespace Foam
{

autoPtr<tcLibraryInterface> tcLibraryInterface::New(const fvMesh& mesh)
{
    IOdictionary dict
    (
        IOobject
        (
            "thermophysicalProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        )
    );

    const word libName =
        dict.lookupOrDefault<word>("library", "neTCLib");

    Info<< "tcLibraryInterface::New : library = " << libName << nl << endl;

    if (libName == "mutationpp")
    {
        FatalErrorInFunction
            << "library 'mutationpp' selected but Mutation++ backend "
            << "is not yet implemented." << nl
            << "Set 'library' to 'neTCLib' or 'default'."
            << exit(FatalError);

        return autoPtr<tcLibraryInterface>(nullptr);  // unreachable
    }
    else if
    (
        libName == "neTCLib"
     || libName == "default"
    )
    {
        // ne2TReactionThermo inherits tcLibraryInterface.
        // Releasing the raw pointer into autoPtr<tcLibraryInterface> is
        // safe because the base destructor is virtual.
        autoPtr<ne2TReactionThermo> native =
            ne2TReactionThermo::New(mesh);

        return autoPtr<tcLibraryInterface>(native.release());
    }
    else
    {
        FatalErrorInFunction
            << "Unknown library value: '" << libName << "'." << nl
            << "Valid options: 'neTCLib', 'default', 'mutationpp'."
            << exit(FatalError);

        return autoPtr<tcLibraryInterface>(nullptr);  // unreachable
    }
}

} // End namespace Foam
