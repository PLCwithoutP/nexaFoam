IOdictionary thermoDict_
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

const dictionary& mix = thermoDict_.subDict("mixture");
const dictionary& sp  = mix.subDict("specie");

word specieKind = Foam::word("");
sp.readIfPresent("type", specieKind);  // SAFE: no fatal if missing

const bool vibrationalCheck = (specieKind == "molecular");
