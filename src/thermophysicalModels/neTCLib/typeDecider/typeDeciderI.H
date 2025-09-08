#ifndef typeDeciderI_H
#define typeDeciderI_H

template<class OneT, class TwoT>
inline bool Foam::typeDecider<OneT,TwoT>::construct
(
    autoPtr<OneT>& p1T,
    autoPtr<TwoT>& p2T,
    const fvMesh& mesh
)
{
    bool check2T = (tempModel == word("heNe2TThermo"));
    bool checkT = (tempModel == word("heNeTThermo"));

    if (check2T)
    {
        p2T = TwoT::New(mesh);
        Info << "Selecting two-temperature model " << nl << endl;
        return true;
    }
    else if (checkT)
    {
        p1T = OneT::New(mesh);
        Info << "Selecting one-temperature model " << nl << endl;
        return false;
    }
    else
    {
        Info
            << "thermoType.type = " << tempModel
            << " should be " << checkT << " or " << check2T
            << endl;
        return false;
    }
}

#endif
