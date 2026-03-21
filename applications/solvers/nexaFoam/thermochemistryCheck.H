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
        else if (mixtureType == "pure2TMixture")
        {
            mixtureCheck = false;
            Info << "Pure mixture is detected." << nl
                 << endl;
        }
        else
        {
            WarningInFunction
                << "Unknown thermoType.mixture = " << mixtureType
                << ". Expected 'multiComponent2TMixture' or 'pure2TMixture'. "
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
