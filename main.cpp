#include "functions.h"

// Declare global variables
const int U = 100;
areaParams aAreaParams;
int iSeed = int(time(0));
default_random_engine generator(iSeed);
cellHolders cell;
epitheliumLocations aEpLocation;

int main()
{
    aAreaParams.iSize = U;
    aAreaParams.dMesenchymeDensity = 0.1;
    aAreaParams.iNumMesenchyme = fNumMesenchyme(aAreaParams.iSize,aAreaParams.dMesenchymeDensity);
    MatrixXi mArea = fCreateInitialArea(aAreaParams.iSize,250,U/2,U/2,aAreaParams.dMesenchymeDensity,cell,aEpLocation);
    MatrixXd mGDNF = fFieldUpdate(mArea, cell, aAreaParams);

    return 0;
}

double fRandUniformDouble()
{
    uniform_real_distribution<double> uDistribution(0.0,1.0);
    double dRand = uDistribution(generator);
    return dRand;
}

