#include "functions.h"

// Declare global variables
extern const int U = 100;
areaParams aAreaParams;
cellHolders cell;
epitheliumLocations aEpLocation;
totals aTotals;
ostringstream os;ofstream GDNF, CellLocations;

int main()
{
    aAreaParams.iSize = U;
    aAreaParams.dMesenchymeDensity = 0.1;
    aAreaParams.iNumMesenchyme = fNumMesenchyme(aAreaParams.iSize,aAreaParams.dMesenchymeDensity);
    MatrixXi mArea = fCreateInitialArea(aAreaParams.iSize,100,U/2,U/2,aAreaParams.dMesenchymeDensity,cell,aEpLocation,aTotals);
    MatrixXd mGDNF = fFieldUpdate(mArea, cell, aAreaParams);
    vector<MatrixXi> aMovieHolderCell;
    vector<MatrixXd> aMovieHolderGDNF;

    int iFrameRate = 2;
    int iT = 100;
    for (int t = 0; t < iT; ++t)
    {
        cout<<t<<"\n";
        fUpdateEpithelium(mGDNF,mArea,cell,aAreaParams,aEpLocation,aTotals);
        if (t%iFrameRate == 0)
        {
            mGDNF = fFieldUpdate(mArea,cell,aAreaParams);
            aMovieHolderCell.push_back(mArea);
            aMovieHolderGDNF.push_back(mGDNF);
        }
    }

    string fileNameCell = "Results-CellLocations.txt";
    string fileNameGDNF = "Results-GDNF.txt";
    WriteVectorToFile(aMovieHolderCell,fileNameCell, iT/iFrameRate);
    WriteVectorToFile(aMovieHolderGDNF,fileNameGDNF, iT/iFrameRate);



    return 0;
}
