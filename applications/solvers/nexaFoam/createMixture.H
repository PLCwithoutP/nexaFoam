autoPtr<ne2TReactionThermo> thermoSelector
(
    ne2TReactionThermo::New(mesh)
);

ne2TReactionThermo& thermo2T = thermoSelector();

auto& mix = thermo2T.composition();

const speciesTable& spec  = mix.species();
