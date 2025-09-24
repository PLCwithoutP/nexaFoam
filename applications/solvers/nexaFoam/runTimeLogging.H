myFile << "Time is : " << runTime.timeName() << endl;
myFile << "Average Translational Cv value in internal cells is : \n" << (thermo2T.CvT())().average().value() << endl;
myFile << "Average Rotational Cv value in internal cells is : \n" << (thermo2T.CvR())().average().value() << endl;
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

Info << "Average rotational internal energy is: " << eR.average().value() << nl;
myFile << "Average rotational internal energy is: " << eR.average().value() << nl;

if (vibrationalCheck)
{
    Info << "Average vibrational internal energy is: " << eV.average().value() << nl;
    myFile << "Average vibrational internal energy is: " << eV.average().value() << nl;
}
