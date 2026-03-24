autoPtr<ne2TReactionThermo> thermoSelector
(
    ne2TReactionThermo::New(mesh)
);

ne2TReactionThermo& thermo2T = thermoSelector();

auto& mix = thermo2T.composition();
const speciesTable& spec = mix.species();

const bool use2T = thermo2T.twoTemperature();
const bool useChemistry = thermo2T.chemistry();
const bool chemistryCapable = thermo2T.chemistryCapable();

const fileName& speciesListFile = thermo2T.speciesListFile();
const fileName& reactionsListFile = thermo2T.reactionsListFile();
const word& chemistryReader = thermo2T.chemistryReader();

Info<< "twoTemperature = " << (use2T ? "true" : "false") << nl
    << "chemistry     = " << (useChemistry ? "true" : "false") << nl
    << "speciesList   = " << speciesListFile << nl;

if (useChemistry)
{
    Info<< "chemistryReader = " << chemistryReader << nl
        << "reactionsList   = " << reactionsListFile << nl;
}

Info<< endl;

if (useChemistry && !chemistryCapable)
{
    FatalErrorInFunction
        << "chemistry=true was requested in thermophysicalProperties, "
        << "but the selected thermo package does not support chemistry."
        << nl
        << "Choose a reaction-capable thermo or set chemistry=false."
        << exit(FatalError);
}

autoPtr<BasicNeChemistryModel<ne2TReactionThermo>> chemistryPtr;

if (useChemistry)
{
    chemistryPtr = BasicNeChemistryModel<ne2TReactionThermo>::New(thermo2T);

    Info<< "Chemistry model created." << nl << endl;
}
