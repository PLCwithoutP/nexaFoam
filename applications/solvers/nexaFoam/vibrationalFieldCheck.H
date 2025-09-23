IOdictionary thermoDict_
(
    IOobject
    (
        "thermophysicalProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    )    
);

const scalar& thetaVib_ = thermoDict_.subDict("mixture").subDict("thermodynamics").get<scalar>("ThetaVibrational");

const bool vibrationalCheck = (thetaVib_ != 1e-6);