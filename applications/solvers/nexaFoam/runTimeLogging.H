myFile << "Time is : " << runTime.timeName() << endl;

forAll(spec, i)
{
    const volScalarField& Yi = Y_species[i];
    const scalar Wi = mix.W(i);
    const scalar mui = mix.mu(i,TTR.average().value());

    const scalar volAvg =
        Yi.weightedAverage(mesh.V()).value();

    Info << "Average mass fraction of " << Yi.name() << " is:  " << volAvg << nl;
    myFile << "Average mass fraction of " << Yi.name() << " is:  " << volAvg << nl;

    Info << "Mol weight of " << spec[i] << " is:  " << Wi << nl;
    myFile << "Mol weight of " << spec[i] << " is:  " << Wi << nl;

    Info << "Average viscosity of " << spec[i] << " is:  " << mui << nl;
    myFile << "Average viscosity of " << spec[i] << " is:  " << mui << nl;
}


myFile << "Average Translational Cv value in internal cells is : \n" << (thermo2T.CvT())().average().value() << endl;
if(vibrationalCheck)
{
    myFile << "Average Rotational Cv value in internal cells is : \n" << (thermo2T.CvR())().average().value() << endl;
}
myFile << "Average Translational-Rotational Cp value in internal cells is : \n" << (thermo2T.CpTR())().average().value() << endl;
if(vibrationalCheck)
{
    myFile << "Average Vibrational Cv value in internal cells is : \n" << (thermo2T.CvVib())().average().value() << endl;
}
myFile << "Average Translational-Rotational mu value in internal cells is : \n" << (thermo2T.mu())().average().value() << endl;
myFile << "Average Translational-Rotational kappa value in internal cells is : \n" << (thermo2T.kappaTR())().average().value() << endl;
if(vibrationalCheck)
{
    myFile << "Average Vibrational kappa value in internal cells is : \n" << (thermo2T.kappaVib())().average().value() << endl;
}

Info << "Average translational-rotational temperature is: " << TTR.average().value() << nl;
myFile << "Average translational-rotational temperature is: " << TTR.average().value() << nl;

if (vibrationalCheck)
{
    Info << "Average vibrational temperature is: " << TVib.average().value() << nl;
    myFile << "Average vibrational temperature is: " << TVib.average().value() << nl;
}

Info << "Average translational-rotational enthalpy is: " << hTR.average().value() << nl;
myFile << "Average translational-rotational enthalpy is: " << hTR.average().value() << nl;

Info << "Average translational internal energy is: " << eT.average().value() << nl;
myFile << "Average translational internal energy is: " << eT.average().value() << nl;

if (vibrationalCheck)
{
    Info << "Average rotational internal energy is: " << eR.average().value() << nl;
    myFile << "Average rotational internal energy is: " << eR.average().value() << nl;
}

if (vibrationalCheck)
{
    Info << "Average vibrational internal energy is: " << eV.average().value() << nl;
    myFile << "Average vibrational internal energy is: " << eV.average().value() << nl;
}
