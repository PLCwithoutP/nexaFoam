
autoPtr<neTThermo>  p1T;
autoPtr<ne2TThermo> p2T;

Foam::typeDecider<neTThermo, ne2TThermo> typeDecideObj(mesh);

bool isIt2Temp = typeDecideObj.construct(p1T, p2T, mesh);
