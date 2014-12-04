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
#include <algorithm>
#include <functional>

using namespace std;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SparseMatrix;
using Eigen::ConjugateGradient;
using Eigen::SimplicialLLT;
using Eigen::SimplicialLDLT;
using Eigen::BiCGSTAB;
using Eigen::SparseLU;


int iSeed = int(time(0));
default_random_engine generator(iSeed);


typedef SparseMatrix<double> SpMatrix;

// Pre-declared
class mesenchyme;
class epithelium;
struct epitheliumLocations;
struct cellHolders;
struct areaParams;
struct totals;
int mod(int,int);
void fCreateRandomEpithelium(int, int, int, MatrixXi &, epitheliumLocations &, cellHolders &,totals &);
void CreateEpithelium(int, int, MatrixXi &, epitheliumLocations &, cellHolders &, totals &);
void CreateMesenchyme(int, int, MatrixXi &, cellHolders &, totals &);
int getConnectivity(int, int, MatrixXi &);
inline int fProbabilitySwitch(double aProbability);
pair<int,MatrixXi> fAllowedEpithelium(const int, const int, const int, MatrixXi &, const areaParams &);
double normcdf(double);
int fMoveOrProlif (int, int, MatrixXi &, const MatrixXd &, const areaParams &);
void fMoveProlifEpithelium(int, int, int,MatrixXi &, const areaParams &, epithelium &, epitheliumLocations &,cellHolders &, totals &);
int RandomInteger(int,int);
int isSamePoint(int, int, int, int);
SpMatrix Laplace(int);
SpMatrix LaplaceNeumann(int);

extern const int U;
static SpMatrix mA = LaplaceNeumann(U);
SimplicialLLT<SparseMatrix<double>> solver;
int iNearestNeighbours = 8;

void solveLinearSystem()
{
    solver.analyzePattern(mA);   // for this step the numerical values of A are not used
}


struct areaParams
{
    int iSize;
    double dMesenchymeDensity;
    int iNumMesenchyme;
    double dDg = 100.0;
    double dGamma = 1;
    double dPMoveVsPProlif = 0.5;
    double dMoveProlifC0 = -10;
    double dMoveProlifC1 = 7;
    int iChemotaxis = 0;
};

struct totals
{
    int iNumEpithelium = 0;
    int iNumMesenchyme = 0;
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
        mesenchyme();
        void Move(int, int, int, int, MatrixXi &);
        int getX();
        int getY();
};

mesenchyme::mesenchyme()
{
    iX = -1;
    iY = -1;
}

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

void mesenchyme::Move(int aX, int aY, int bX, int bY, MatrixXi &aArea)
{
    //The epithelium move function should already have taken care of that square, so no need to set it to zero
    this->iX=bX;
    this->iY=bY;
    aArea(bX,bY) = -1;
}

class epithelium
{
    private:
        int iX, iY;

    public:

        epithelium(int,int,MatrixXi &, epitheliumLocations &);
        int getX();
        int getY();
        void Move(int, int, int, int, MatrixXi &, totals &);
        void Prolif(int, int, int, int, MatrixXi &,epitheliumLocations &, cellHolders &, totals &);
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


void epithelium::Move(int aX, int aY, int bX, int bY, MatrixXi &aArea, totals &aTotals)
{
    this->iX=bX;
    this->iY=bY;
    aArea(aX,aY) = 0;
    aArea(bX,bY) = 1;
}

void epithelium::Prolif(int aX, int aY, int bX, int bY, MatrixXi &aArea, epitheliumLocations &aEpitheliumLocation, cellHolders &aCell, totals &aTotals)
{
    aArea(bX,bY) = 1;
    CreateEpithelium(bX, bY, aArea, aEpitheliumLocation, aCell,aTotals);
}

// Functions
double fRandUniformDouble();

int fNumMesenchyme(int iSize, double dMesenchymeDensity)
{
    int iNumMesenchyme;
    iNumMesenchyme = int(dMesenchymeDensity*pow(iSize,2));
    return iNumMesenchyme;
}


MatrixXi fCreateInitialArea(const int iSize, const int iNumEpithelium, const int iXstart, const int iYstart, const double dMesenchymeDensity, cellHolders &aCell, epitheliumLocations &aEpitheliumLocation, totals &aTotals)
{
    solveLinearSystem();
    int iNumMesenchyme = fNumMesenchyme(iSize,dMesenchymeDensity);
    MatrixXi aArea = MatrixXi::Zero(iSize,iSize);
    fCreateRandomEpithelium(iNumEpithelium, iXstart, iYstart, aArea, aEpitheliumLocation, aCell, aTotals);

    //Pick random locations for the initial mesenchyme cells
    double iRandX, iRandY;
    int i = 0;
    while(i < iNumMesenchyme)
    {
        iRandX = int(fRandUniformDouble()*iSize);
        iRandY = int(fRandUniformDouble()*iSize);
        if (aArea(iRandX,iRandY)==0)
        {
            CreateMesenchyme(iRandX, iRandY, aArea, aCell, aTotals);
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
        for (int i = 0; i < 4; i++)
        {   //Rather than make a copy, we will check to see that the connected cell is not the original (as it has moved)
            if (mArea(mIndices(i,0),mIndices(i,1))== 1 && (mIndices(i,0)!=aX && mIndices(i,1)!=aY))
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
    MatrixXi mIndices;
    if (iNearestNeighbours==4)
    {
        mIndices = fAllIndices8Neighbours(aX,aY,aArea);
    }
    else
    {
        mIndices = fAllIndices8Neighbours(aX,aY,aArea);
    }

    MatrixXi mAllowed = MatrixXi::Zero(iNearestNeighbours,2);

    int k = 0;
    for (int i = 0; i < iNearestNeighbours; i++)
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
        if (iMove == 0)
        {

            if (fActiveConnected(aX, aY, mAllowed(i,0), mAllowed(i,1), 0, aArea)==1)
            {
                mFeasible(k,0) = mAllowed(i,0);
                mFeasible(k,1) = mAllowed(i,1);
                k++;

            }
        }
        else
        {

            if (fActiveConnected(aX, aY, mAllowed(i,0), mAllowed(i,1), 1, aArea)==1)
            {

                mFeasible(k,0) = mAllowed(i,0);
                mFeasible(k,1) = mAllowed(i,1);
                k++;

            }
        }

    }

return mAllowed;

}

void fCreateRandomEpithelium(const int iNumEpithelium, int iXstart, int iYstart, MatrixXi &aArea, epitheliumLocations &aEpitheliumLocation, cellHolders &aCell, totals &aTotals)
{
    // Create an epithelium cell in the middle
    CreateEpithelium(iXstart,iYstart,aArea,aEpitheliumLocation,aCell,aTotals);
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
            CreateEpithelium(iXstart,iYstart,aArea,aEpitheliumLocation,aCell,aTotals);
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
                    CreateEpithelium(mIndices(iRand,0),mIndices(iRand,1),aArea,aEpitheliumLocation,aCell,aTotals);
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
void CreateEpithelium(const int aX, const int aY, MatrixXi &aArea, epitheliumLocations &aEpitheliumLocation, cellHolders &aCell, totals &aTotals)
{
    epithelium aEpithelium(aX,aY,aArea,aEpitheliumLocation);
    aTotals.iNumEpithelium++;
    aCell.groupEpithelium.push_back(aEpithelium);
}

// A function which creates an epithelium and adds it to the list of epithelium
void CreateMesenchyme(const int aX, const int aY, MatrixXi &aArea, cellHolders &aCell, totals &aTotals)
{
    mesenchyme aMesenchyme(aX,aY,aArea);
    aTotals.iNumMesenchyme++;
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

MatrixXd fFieldUpdate(const MatrixXi &aArea, const cellHolders &aCell, const areaParams & aAreaParams)
{

    int iSize = aArea.rows();
    double dDg = aAreaParams.dDg;
    double dGamma = aAreaParams.dGamma;
    SpMatrix A = mA;

    // Go through the epithelium and amend the Laplacian in the correct spots
    vector<epithelium> aVector = aCell.groupEpithelium;
    int iX, iY;
    for (std::vector<epithelium>::iterator it = aVector.begin() ; it != aVector.end(); ++it)
    {
        epithelium aEpithelium = *it;
        iX = aEpithelium.getX();
        iY = aEpithelium.getY();
        A.coeffRef(iX*iSize+iY, iX*iSize+iY) += (1.0/dDg);
    }

    // Now create the rhs 'b' vector and modify it according to mesenchyme locations
    VectorXd b = VectorXd::Zero(iSize*iSize);
    vector<mesenchyme> bVector = aCell.groupMesenchyme;
    for (std::vector<mesenchyme>::iterator it = bVector.begin() ; it != bVector.end(); ++it)
    {
        mesenchyme aMesenchyme = *it;
        iX = aMesenchyme.getX();
        iY = aMesenchyme.getY();
        b(iX*iSize + iY) = (dGamma/dDg);
    }

    // Now solving for the GDNF concentration
        solver.factorize(A);
        MatrixXd vGDNF = solver.solve(b);
//        solver.compute(A).solve(b);
//        MatrixXd vGDNF = solver.compute(A).solve(b);


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

void fUpdateEpithelium(const MatrixXd &mGDNF, MatrixXi &aArea, cellHolders &aCell, const areaParams & aAreaParams,epitheliumLocations &aEpitheliumLocations,totals &aTotals)
{
    double dPMoveVsPProlif = aAreaParams.dPMoveVsPProlif;
    double dMoveProlifC0 = aAreaParams.dMoveProlifC0;
    double dMoveProlifC1 = aAreaParams.dMoveProlifC1;
    int iChemotaxis = aAreaParams.iChemotaxis;

    //Create a randomly sorted vector of the epithelium
    vector<epithelium> & aVector = aCell.groupEpithelium;
    random_shuffle(aVector.begin(),aVector.end());
    int iX, iY, iMove;
    int iLength = aVector.size();
    //Use a static for loop rather than an iterator, since we are changing its length
    for (int i = 0; i < iLength; ++i)
    {
        epithelium & aEpithelium = aVector[i];
        iX = aEpithelium.getX();
        iY = aEpithelium.getY();

        //Select between a move or prolif probabilistically
        iMove = fProbabilitySwitch(dPMoveVsPProlif);

        //Determine whether or not there are allowed moves for the epithelium, taking into account the mesenchyme
        pair<int,MatrixXi> aPairAllowed = fAllowedEpithelium(iX,iY,iMove,aArea,aAreaParams);

        //Allow a move/prolif to depend on the local concentration of GDNF
        if (aPairAllowed.first>0)
        {
            int cEvent = fMoveOrProlif(iX,iY,aArea,mGDNF,aAreaParams);

            //Select between the moves either randomly or dependent on the local concentration of GDNF
            if (cEvent>0)
            {
                fMoveProlifEpithelium(iX,iY,iMove,aArea,aAreaParams,aEpithelium,aEpitheliumLocations,aCell,aTotals);
            }
        }
    }
}

// Returns a matrix of indices for those cells around the cell location which are feasible
// meaning either they are vacant/mesenchyme and active
MatrixXi fAllowedActiveEpithelium(const int aX, const int aY, const int iMove, MatrixXi &aArea)
{
    MatrixXi mAllowed;
    if (iNearestNeighbours == 4)
    {
        mAllowed = fAllIndices4Neighbours(aX,aY,aArea);
    }
    else
    {
        mAllowed = fAllIndices8Neighbours(aX,aY,aArea);
    }

    MatrixXi mFeasible = MatrixXi::Zero(iNearestNeighbours,2);
    int k = 0;

    for (int i = 0; i < iNearestNeighbours; ++i)
    {
        if (iMove == 1)
        {
            if (fActiveConnected(aX, aY, mAllowed(i,0), mAllowed(i,1), 1, aArea)==1)
            {
                mFeasible(k,0) = mAllowed(i,0);
                mFeasible(k,1) = mAllowed(i,1);
                k++;
            }
        }
        else
        {
            if (fActiveConnected(aX, aY, mAllowed(i,0), mAllowed(i,1), 0, aArea)==1)
            {
                mFeasible(k,0) = mAllowed(i,0);
                mFeasible(k,1) = mAllowed(i,1);
                k++;
            }
        }

    }
return mAllowed;
}


//Determines whether moves/prolifs are allowed, taking into account the mesenchyme.
pair<int,MatrixXi> fAllowedEpithelium(const int iX, const int iY, const int iMove, MatrixXi &aArea, const areaParams & aAreaParams)
{
    pair<int,MatrixXi> aPair;
    int cAllowed = 0;
    MatrixXi mAllIndices = fAllowedActiveEpithelium(iX, iY, iMove, aArea);
    int iRows = mAllIndices.rows();
    MatrixXi mAllowedIndices = MatrixXi::Zero(iRows,2);

    int k = 0;
    for (int i = 0; i < iRows; ++i)
    {
        if (aArea(mAllIndices(i,0),mAllIndices(i,1)) == 0)//Vacant
            {
                mAllowedIndices(k,0) = mAllIndices(i,0);
                mAllowedIndices(k,1) = mAllIndices(i,1);
                k++;
            }
        else if (aArea(mAllIndices(i,0),mAllIndices(i,1)) == -1)//Mesenchyme
        {
            MatrixXi mTemp = fFindAnyVacant(mAllIndices(i,0),mAllIndices(i,1),aArea);
            if (mTemp.rows()>0)
            {
                mAllowedIndices(k,0) = mAllIndices(i,0);
                mAllowedIndices(k,1) = mAllIndices(i,1);
                k++;
            }
        }
    }
    cAllowed = k;
    aPair = make_pair(cAllowed,mAllowedIndices);
    return aPair;
}

double fRandUniformDouble()
{
    uniform_real_distribution<double> uDistribution(0.0,1.0);
    double dRand = uDistribution(generator);
    return dRand;
}

// Determines whether an event takes place given its prior probability
inline int fProbabilitySwitch(double aProbability)
{
    int cEvent;
    double aRand = fRandUniformDouble();
    if (aProbability > aRand) cEvent = 1;
    else cEvent = 0;
    return cEvent;
}

// Determines whether a move or proliferation takes place, dependent on the local probability of GDNF
int fMoveOrProlif (int iX, int iY, MatrixXi &aArea, const MatrixXd &mGDNF, const areaParams & aAreaParams)
{
    // Calculate the 'X' coordinate of the normal
    double dXNormal = aAreaParams.dMoveProlifC0 + aAreaParams.dMoveProlifC1*mGDNF(iX,iY);

    // Calculate probability from normal CDF
    double dXNormalProb = normcdf(dXNormal);

    // Test it
    int cEvent = fProbabilitySwitch(dXNormalProb);
    return cEvent;
}

double normcdf(double x)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

struct IsSameMesenchyme: public std::binary_function< mesenchyme, pair<int,int>, bool > {
  bool operator () (mesenchyme &aMesenchyme, const pair<int,int> &pairLocation) const {
    int aX = aMesenchyme.getX(); int aY = aMesenchyme.getY();
    int bX = pairLocation.first; int bY = pairLocation.second;
    bool aBool = (aX==bX && aY==bY);
//    if (aBool)
//        cout<<aX<<"  "<<aY<<"  "<<bX<<"   "<<bY<<"\n";
    return aBool;
    }
  };

void fMoveProlifEpithelium(int iX, int iY, int iMove,MatrixXi &aArea, const areaParams &aAreaParams, epithelium &aEpithelium, epitheliumLocations &aEpitheliumLocations,cellHolders &aCell,totals &aTotals)
{
    pair<int,MatrixXi> aPairAllowed = fAllowedEpithelium(iX,iY,iMove,aArea,aAreaParams);
    if (aPairAllowed.first == 0) {cout<<"An error has been made in fMoveProlifEpithelium"<<"\n";return;};
    //Choose a position randomly
    int aRand, bX, bY;
    aRand = RandomInteger(0,aPairAllowed.first);
    MatrixXi mTemp = aPairAllowed.second;
    bX = mTemp(aRand,0);
    bY = mTemp(aRand,1);
    int iCellType = aArea(bX,bY);
    if (iCellType==1){cout<<"An error has been made where an epithelium is being moved into in fMoveProlifEpithelium"<<"\n";return;};

    //Moving
    if (iMove == 1)
    {
        if (iCellType == 0) //Vacant spot
        {
            aEpithelium.Move(iX,iY,bX,bY,aArea,aTotals);
            return;
        }
        else if (iCellType == -1)
        {
            //Proliferate the epithelium and choose a vacant space for the mesenchyme to move into
            aEpithelium.Move(iX,iY,bX,bY,aArea,aTotals);
            MatrixXi mMesenchymeVacant = fFindAnyVacant(bX,bY,aArea);
            if (mMesenchymeVacant.rows()==0){cout<<"Something has gone wrong in fMoveProlifEpithelium where a mesenchyme has nowhere to go\n";};
            int aRandInt = RandomInteger(0,mMesenchymeVacant.rows());

            //Find the relevant mesenchyme by searching for it in a list of mesenchyme. Have used a binary adaptable predicate
            //by creating a struct which returns True if all location parameters match).
            //This is to get round the fact that the predicate of find_if is by its nature unary.
            //A predicate is a function returning a boolean.
            vector<mesenchyme> & aMesenchymeGroup = aCell.groupMesenchyme;
            pair<int,int> aPair = make_pair(bX,bY);
            vector<mesenchyme>::iterator it = find_if(aMesenchymeGroup.begin(),aMesenchymeGroup.end(), bind2nd(IsSameMesenchyme(),aPair));
            mesenchyme aMesenchyme = *it;
            aMesenchyme.Move(bX,bY,mMesenchymeVacant(aRandInt,0),mMesenchymeVacant(aRandInt,1),aArea);
            return;
        }
    }
    //Proliferating into a vacant spot
    else
    {
        if (iCellType == 0)
        {
            aEpithelium.Prolif(iX,iY,bX,bY,aArea,aEpitheliumLocations,aCell,aTotals);
            return;
        }
        else if (iCellType == -1)
        {
            //Proliferate the epithelium and choose a vacant space for the mesenchyme to move into
            aEpithelium.Prolif(iX,iY,bX,bY,aArea,aEpitheliumLocations,aCell,aTotals);
            MatrixXi mMesenchymeVacant = fFindAnyVacant(bX,bY,aArea);
            if (mMesenchymeVacant.rows()==0){cout<<"Something has gone wrong in fMoveProlifEpithelium where a mesenchyme has nowhere to go\n";};
            int aRandInt = RandomInteger(0,mMesenchymeVacant.rows());

            //Find the relevant mesenchyme by searching for it in a list of mesenchyme. Have used a binary adaptable predicate
            //by creating a struct which returns True if all location parameters match).
            //This is to get round the fact that the predicate of find_if is by its nature unary.
            //A predicate is a function returning a boolean.
            vector<mesenchyme> & aMesenchymeGroup = aCell.groupMesenchyme;
            pair<int,int> aPair = make_pair(bX,bY);
            vector<mesenchyme>::iterator it = find_if(aMesenchymeGroup.begin(),aMesenchymeGroup.end(), bind2nd(IsSameMesenchyme(),aPair));
            mesenchyme aMesenchyme = *it;
            aMesenchyme.Move(bX,bY,mMesenchymeVacant(aRandInt,0),mMesenchymeVacant(aRandInt,1),aArea);
            return;
        }
    }
}


int RandomInteger(int iMin, int iMax)
{
    int aRandInt = int(fRandUniformDouble()*(iMax-iMin) + iMin);
    return aRandInt;
}

template<typename T>
void WriteVectorToFile(vector<T> aVector ,string aFilename, int iNumComponents)
{
    ostringstream os;
    ofstream fileTemp;
    os<<aFilename;
    fileTemp.open(os.str().c_str());
    for (int t = 0; t < iNumComponents; ++t)
    {
        fileTemp<<aVector[t]<<"\n";
    }
    fileTemp.close();os.str("");
}

template<typename T>
void WriteToFile(T aThing ,string aFilename)
{
    ostringstream os;
    ofstream fileTemp;
    os<<aFilename;
    fileTemp.open(os.str().c_str());
    fileTemp<<aThing;
    fileTemp.close();os.str("");
}

SpMatrix LaplaceNeumann(int m)
{
    const int n = m*m, n0 = n - 1, n_m = n - m, m0 = m - 1,
    NumNonZeros = 5*n;

    SpMatrix A(n,n);
    A.reserve(Eigen::VectorXi::Constant(n,5));
//    cout<<NumNonZeros<<"\n";

    for (int j = 0; j < n; j++)
    {
        //Diagonals
        A.insert(j,j) = 4;
        if (j > 0 && (j < m0 || j > n - m)) //First block diagonal
        {
            A.coeffRef(j,j) = 3;
        }
        if (j == m0 || j == n - m)
        {
            A.coeffRef(j,j) = 2;
        }
        if (j > m0 && (j+1)%m == 0)
        {
            A.coeffRef(j,j) = 3;
        }
        if (j > m0 && j%m == 0 && j < n - m)
        {
            A.coeffRef(j,j) = 3;
        }

        if (j == 0 || j == n0) //First and last
        {
            A.coeffRef(j,j) = 2;
        }

        //Near diagonals
        if (j > 0 && j%m!=0)
        {
            A.insert(j,j-1) = -1;
            A.insert(j-1,j) = -1;
        }

        //Far diagonals
        if (j >= m)
        {
            A.insert(j,j-m) = -1;
            A.insert(j-m,j) = -1;
        }

//        if (j >= m)
//            {
//                A.insert(j,j-m) = c; A.insert(j-m,j) = c;
//            }
//        if (j > 0 && j%m != 0) A.insert(j-1,j) = b;
//        A.insert(j,j) = d;
//        if (j < n0 && (j+1)%m != 0) A.insert(j+1,j) = b;

    }

return A;
}

//MatrixXd fSolveSystemOld()

#endif // FUNCTIONS_H_INCLUDED
