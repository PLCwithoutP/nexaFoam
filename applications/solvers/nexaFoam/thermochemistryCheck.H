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

// --- Two-temperature check
bool vibrationalCheck = false;

if (thermoDict_.found("twoTemperature"))
{
    const word twoTemperatureValue(thermoDict_.lookup("twoTemperature"));

    if (twoTemperatureValue == "true")
    {
        vibrationalCheck = true;

        Info<< "twoTemperature = true. Vibrational mode is enabled."
            << nl << endl;
    }
    else if (twoTemperatureValue == "false")
    {
        vibrationalCheck = false;

        Info<< "twoTemperature = false. Vibrational mode is disabled."
            << nl << endl;
    }
    else
    {
        FatalErrorInFunction
            << "Invalid value for 'twoTemperature' in thermophysicalProperties: "
            << twoTemperatureValue << nl
            << "Expected 'true' or 'false'."
            << exit(FatalError);
    }
}
else
{
    FatalErrorInFunction
        << "Missing required entry 'twoTemperature' in thermophysicalProperties."
        << nl
        << "Add either:" << nl
        << "    twoTemperature true;" << nl
        << "or" << nl
        << "    twoTemperature false;"
        << exit(FatalError);
}

// --- Mixture-type check (multiComponent2TMixture vs pureMixture)
bool mixtureCheck = false;

if (thermoDict_.found("thermoType"))
{
    const dictionary& thermoTypeDict = thermoDict_.subDict("thermoType");

    if (thermoTypeDict.found("mixture"))
    {
        const word mixtureType(thermoTypeDict.lookup("mixture"));

        if (mixtureType == "multiComponent2TMixture")
        {
            mixtureCheck = true;
            Info << "Multiple species are found in mixture." << nl
                 << endl;
        }
        else if (mixtureType == "pureMixture")
        {
            mixtureCheck = false;
            Info << "Pure mixture is detected." << nl
                 << endl;
        }
        else
        {
            WarningInFunction
                << "Unknown thermoType.mixture = " << mixtureType
                << ". Expected 'multiComponent2TMixture' or 'pureMixture'. "
                << "Defaulting mixtureCheck=false." << nl << endl;

            mixtureCheck = false;
        }
    }
    else
    {
        WarningInFunction
            << "thermoType dictionary found, but 'mixture' entry missing. "
            << "Defaulting mixtureCheck=false." << nl << endl;
    }
}
else
{
    WarningInFunction
        << "'thermoType' dictionary missing in thermophysicalProperties. "
        << "Defaulting mixtureCheck=false." << nl << endl;
}
