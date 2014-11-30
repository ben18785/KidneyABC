#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <random>
#include <time.h>
#include <list>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/SparseExtra>


using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;
using Eigen::SimplicialLLT;
using Eigen::SimplicialLDLT;
using Eigen::BiCGSTAB;
using Eigen::SparseLU;

using namespace std;

typedef SparseMatrix<double> SpMatrix;

// Pre-declared
class mesenchyme;
class epithelium;
struct epitheliumLocations;
struct cellHolders;
int mod(int,int);
void fCreateRandomEpithelium(int, int, int, MatrixXi &, epitheliumLocations &, cellHolders &);
void CreateEpithelium(int, int, MatrixXi &, epitheliumLocations &, cellHolders &);
void CreateMesenchyme(int, int, MatrixXi &, cellHolders &);
int getConnectivity(int, int, MatrixXi &);


struct areaParams
{
    int iSize;
    double dMesenchymeDensity;
    int iNumMesenchyme;
    double dDg = 100.0;
    double dGamma = 1;
};

struct epitheliumLocations
{
    vector<MatrixXi> locations;
};


// Holds the cell types
struct cellHolders
{
    vector<mesenchyme> groupMesenchyme;
    vector<epithelium> groupEpithelium;
};

//Classes
class mesenchyme
{
    private:
    int iX, iY;

    public:
        mesenchyme(int,int,MatrixXi &);
        int getX();
        int getY();
};

mesenchyme::mesenchyme(int aX, int aY, MatrixXi &mArea)
{
    iX = aX;
    iY = aY;

    if (mArea(aX,aY)==0 || mArea(aX,aY)==1)
    {
        mArea(aX,aY) = -1;
    }
    else
    {
            cout<<"Tried to create a mesenchyme in a location where there is already a mesenchyme\n";
    }

}

int mesenchyme::getX()
{
    return iX;
}

int mesenchyme::getY()
{
    return iY;
}

class epithelium
{
    private:
    int iX, iY;

    public:
        epithelium(int,int,MatrixXi &, epitheliumLocations &);
        int getX();
        int getY();
};

epithelium::epithelium(int aX, int aY, MatrixXi &mArea, epitheliumLocations &aEpitheliumLocation)
{
    iX = aX;
    iY = aY;

    if (mArea(aX,aY)==0 || mArea(aX,aY)==1)
    {
        mArea(aX,aY) = 1;
        MatrixXi aRowLocation(1,2);
        aRowLocation << aX, aY;
        aEpitheliumLocation.locations.push_back(aRowLocation);

    }
    else
    {
            cout<<"Tried to create an epithelium in a location where there is already a mesenchyme\n";
    }

}

int epithelium::getX()
{
    return iX;
}

int epithelium::getY()
{
    return iY;
}


// Functions
double fRandUniformDouble();

int fNumMesenchyme(int iSize, double dMesenchymeDensity)
{
    int iNumMesenchyme;
    iNumMesenchyme = int(dMesenchymeDensity*pow(iSize,2));
    return iNumMesenchyme;
}


MatrixXi fCreateInitialArea(const int iSize, const int iNumEpithelium, const int iXstart, const int iYstart, const double dMesenchymeDensity, cellHolders &aCell, epitheliumLocations &aEpitheliumLocation)
{
    int iNumMesenchyme = fNumMesenchyme(iSize,dMesenchymeDensity);
    MatrixXi aArea = MatrixXi::Zero(iSize,iSize);
    fCreateRandomEpithelium(iNumEpithelium, iXstart, iYstart, aArea, aEpitheliumLocation, aCell);

    //Pick random locations for the initial mesenchyme cells
    double iRandX, iRandY;
    int i = 0;
    while(i < iNumMesenchyme)
    {
        iRandX = int(fRandUniformDouble()*iSize);
        iRandY = int(fRandUniformDouble()*iSize);
        if (aArea(iRandX,iRandY)==0)
        {
            CreateMesenchyme(iRandX, iRandY, aArea, aCell);
            i++;
        }
    }

    return aArea;
}


MatrixXi fAllIndices8Neighbours(const int aX, const int aY,MatrixXi &aArea)
{
    // Get size
    int iDepth, iWidth;
    iDepth = aArea.rows();
    iWidth = aArea.cols();

    // Initialise the matrix
    MatrixXi mIndices = MatrixXi::Zero(8,2);

    // Now get the neighbours allowing toroidal domain
    // Above
    mIndices(0,0) = mod(aX+1,iDepth);
    mIndices(0,1) = aY;
    // Up right
    mIndices(1,0) = mod(aX+1,iDepth);
    mIndices(1,1) = mod(aY+1,iWidth);
    // Right
    mIndices(2,0) = aX;
    mIndices(2,1) = mod(aY+1,iWidth);
    // Down right
    mIndices(3,0) = mod(aX-1,iDepth);
    mIndices(3,1) = mod(aY+1,iWidth);
    // Down
    mIndices(4,0) = mod(aX-1,iDepth);
    mIndices(4,1) = aY;
    // Down left
    mIndices(5,0) = mod(aX-1,iDepth);
    mIndices(5,1) = mod(aY-1,iWidth);
    // Left
    mIndices(6,0) =  aX;
    mIndices(6,1) = mod(aY-1,iWidth);
    // Left up
    mIndices(7,0) = mod(aX+1,iDepth);
    mIndices(7,1) = mod(aY-1,iWidth);


    return mIndices;
}

MatrixXi fAllIndices4Neighbours(const int aX, const int aY,MatrixXi &aArea)
{

    // Get size
    int iDepth, iWidth;
    iDepth = aArea.rows();
    iWidth = aArea.cols();


    // Initialise the matrix
    MatrixXi mIndices = MatrixXi::Zero(4,2);

    // Now get the neighbours allowing toroidal domain
    // Above

    mIndices(0,0) = mod(aX+1,iDepth);
    mIndices(0,1) = aY;
    // Right
    mIndices(1,0) = aX;
    mIndices(1,1) = mod(aY+1,iWidth);
    // Down
    mIndices(2,0) = mod(aX-1,iDepth);
    mIndices(2,1) = aY;
    // Left
    mIndices(3,0) =  aX;
    mIndices(3,1) = mod(aY-1,iWidth);

    return mIndices;
}



inline int mod(const int iA, const int iB)
{
    return iA-iB*floor(double(iA)/double(iB));
}

// Determines whether the cell position selected is connected. Current position (aX,aY) either
// move or prolif to (bX,bY)
int fActiveConnected(const int aX, const int aY, const int bX, const int bY, const int iMove, MatrixXi &mArea)
{

    MatrixXi mIndices = fAllIndices4Neighbours(bX,bY,mArea);
    int cConnected = 0;

    if (iMove == 1) // Move
    {
        MatrixXi mAreaCopyTemp = mArea;
        mAreaCopyTemp(aX,aY) = 0;
        for (int i = 0; i < 4; i++)
        {

            if (mAreaCopyTemp(mIndices(i,0),mIndices(i,1))== 1)
            {

                cConnected = 1;
                return cConnected;
            }
        }
    }
    else // Prolif
    {
        for (int i = 0; i < 4; i++)
        {

            if (mArea(mIndices(i,0),mIndices(i,1))== 1)
            {

                cConnected = 1;
                return cConnected;
            }
        }

    }

    return cConnected;
}

//Find any vacant cells around a given position
MatrixXi fFindAnyVacant(const int aX, const int aY, MatrixXi &aArea)
{
    MatrixXi mIndices = fAllIndices4Neighbours(aX,aY,aArea);
    MatrixXi mAllowed = MatrixXi::Zero(4,2);

    int k = 0;
    for (int i = 0; i < 4; i++)
    {
        if(aArea(mIndices(i,0),mIndices(i,1))==0)
        {
            mAllowed(k,0) = mIndices(i,0);
            mAllowed(k,1) = mIndices(i,1);
            k++;
        }
    }
    mAllowed.conservativeResize(k,2);
    return mAllowed;
}

// Returns a matrix of indices for those cells around the cell location which are feasible
// meaning both vacant and active
MatrixXi fAllowedActive(const int aX, const int aY, const int iMove, MatrixXi &aArea)
{

    MatrixXi mAllowed = fFindAnyVacant(aX,aY,aArea);


    int iLength = mAllowed.rows();
    MatrixXi mFeasible = MatrixXi::Zero(iLength,2);
    int k = 0;

    if (iLength == 0)
    {
        return mFeasible;
    }

    for (int i = 0; i < iLength; i++)
    {

        if (fActiveConnected(aX, aY, mAllowed(i,0), mAllowed(i,1), 0, aArea)==1)
        {

            mFeasible(k,0) = mAllowed(i,0);
            mFeasible(k,1) = mAllowed(i,1);
            k++;

        }

    }

return mAllowed;

}

void fCreateRandomEpithelium(const int iNumEpithelium, int iXstart, int iYstart, MatrixXi &aArea, epitheliumLocations &aEpitheliumLocation, cellHolders &aCell)
{
    // Create an epithelium cell in the middle
    CreateEpithelium(iXstart,iYstart,aArea,aEpitheliumLocation,aCell);

    MatrixXi mAllowed;
    int cCount = 0;
    int iOldXstart = iXstart;
    int iOldYstart = iYstart;

    int iInitialNumberEpithelium = int(double(iNumEpithelium)*0.5);
    int iRemainingEpithelium = iNumEpithelium - iInitialNumberEpithelium;

    // Create initial (small) clump
    while (cCount < iInitialNumberEpithelium-1)
    {
        mAllowed = fAllowedActive(iXstart,iYstart,0,aArea);
        int iNumRows = mAllowed.rows();
        if (iNumRows>0)
        {
            int iRandRow = int(fRandUniformDouble()*iNumRows);
            int iOldXstart = iXstart;
            int iOldYstart = iYstart;
            iXstart = mAllowed(iRandRow,0);
            iYstart = mAllowed(iRandRow,1);
            CreateEpithelium(iXstart,iYstart,aArea,aEpitheliumLocation,aCell);
            cCount++;
        }
        else
        {
            vector<MatrixXi> aVector = aEpitheliumLocation.locations;
            random_shuffle(aVector.begin(),aVector.end());
            iXstart = aVector[0](0,0);
            iYstart = aVector[0](0,1);
        }

    }

    int iCount = 0;
    while (iCount < iRemainingEpithelium)
    {
        for (int j = 0; j < aEpitheliumLocation.locations.size(); j++)
        {
//            cout<<"j = "<<j<<"\n";
            vector<MatrixXi> aVector = aEpitheliumLocation.locations;
            random_shuffle(aVector.begin(),aVector.end());
            if (getConnectivity(aVector[j](0,0),aVector[j](0,1),aArea) > 1)
            {
                MatrixXi mIndices = fAllowedActive(aVector[j](0,0),aVector[j](0,1),0,aArea);
                if (mIndices.rows()>0)
                {
                    int iRand = int(mIndices.rows()*fRandUniformDouble());
                    CreateEpithelium(mIndices(iRand,0),mIndices(iRand,1),aArea,aEpitheliumLocation,aCell);
                    iCount++;
//                    cout<<iCount<<"\n";
                    if (iCount>iRemainingEpithelium-1)
                    {
                        break;
                    }

                }
            }
        }
    }

}

// Gets the number of connections to other epithelium for a given cell location
int getConnectivity(const int aX, const int aY, MatrixXi &aArea)
{
    MatrixXi mIndices = fAllIndices4Neighbours(aX,aY,aArea);
    int cConnectivityCount = 0;
    for (int i = 0; i < 4; i++)
    {
        if(aArea(mIndices(i,0),mIndices(i,1))==1)
        {
            cConnectivityCount++;
        }
    }
    return cConnectivityCount;
}

// A function which creates an epithelium and adds it to the list of epithelium
void CreateEpithelium(const int aX, const int aY, MatrixXi &aArea, epitheliumLocations &aEpitheliumLocation, cellHolders &aCell)
{
    epithelium aEpithelium(aX,aY,aArea,aEpitheliumLocation);
    aCell.groupEpithelium.push_back(aEpithelium);
}

// A function which creates an epithelium and adds it to the list of epithelium
void CreateMesenchyme(const int aX, const int aY, MatrixXi &aArea, cellHolders &aCell)
{
    mesenchyme aMesenchyme(aX,aY,aArea);
    aCell.groupMesenchyme.push_back(aMesenchyme);
}


SpMatrix Laplace(int m)
{
    const int n = m*m, n0 = n - 1, n_m = n - m,
    NumNonZeros = 5*n;

    SpMatrix A(n,n);
    A.reserve(Eigen::VectorXi::Constant(n,5));
//    cout<<NumNonZeros<<"\n";

    double d = 4, b = -1, c = -1;

    for (int j = 0; j < n; j++)
    {
        if (j >= m)
            {
                A.insert(j,j-m) = c; A.insert(j-m,j) = c;
            }
        if (j > 0 && j%m != 0) A.insert(j-1,j) = b;
        A.insert(j,j) = d;
        if (j < n0 && (j+1)%m != 0) A.insert(j+1,j) = b;

    }

return A;
}

MatrixXd fFieldUpdate(MatrixXi &aArea, cellHolders &aCell, areaParams & aAreaParams)
{
    int iSize = aArea.rows();
    SpMatrix A = Laplace(iSize);

    double dDg = aAreaParams.dDg;
    double dGamma = aAreaParams.dGamma;


    // Go through the epithelium and amend the Laplacian in the correct spots
    vector<epithelium> aVector = aCell.groupEpithelium;
    for (std::vector<epithelium>::iterator it = aVector.begin() ; it != aVector.end(); ++it)
    {
        epithelium aEpithelium = *it;
        int iX = aEpithelium.getX();
        int iY = aEpithelium.getY();
        A.coeffRef(iX*iSize+iY, iX*iSize+iY) += (1.0/dDg);
    }

    // Now create the rhs 'b' vector and modify it according to mesenchyme locations
    VectorXd b = VectorXd::Zero(iSize*iSize);
    vector<mesenchyme> bVector = aCell.groupMesenchyme;
    for (std::vector<mesenchyme>::iterator it = bVector.begin() ; it != bVector.end(); ++it)
    {
        mesenchyme aMesenchyme = *it;
        int iX = aMesenchyme.getX();
        int iY = aMesenchyme.getY();
        b(iX*iSize + iY) = (dGamma/dDg);
    }

    // Now solving for the GDNF concentration
    SimplicialLLT<SparseMatrix<double>> solver;
    MatrixXd vGDNF = solver.compute(A).solve(b);

    // Now going through and folding the vector into a matrix
    MatrixXd mGDNF = MatrixXd::Zero(iSize,iSize);
    for (int i = 0; i < iSize; ++i)
    {
        for (int j = 0; j < iSize; ++j)
        {
            mGDNF(i,j) = vGDNF(i*iSize + j);
        }
    }
    return mGDNF;
}

#endif // FUNCTIONS_H_INCLUDED
