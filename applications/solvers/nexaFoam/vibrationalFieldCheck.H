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

const bool vibrationalCheck = true;
