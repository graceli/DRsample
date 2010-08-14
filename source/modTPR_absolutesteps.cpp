#include <byteswap.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits.h>


static const int nExpectedMagicNumber = 14;

/***************************************************************************//*!
*//****************************************************************************/
void
updateTPR(char* sFile,  int nAdditionalStepLoc, int nAdditionalSteps, 
    int nPullValueLoc, float fPullValue, int nPullKValueLoc,
    float fPullKValue, int nInitStepLoc)
{
    bool bSwap = false;
    int nVal;
    std::fstream tpr(sFile, std::ios::in | std::ios::out | std::ios::binary);

    if (!tpr.is_open())
    {
        std::cerr << "Error opening " << sFile << ". Aborting" << std::endl;
        exit(-1);
    }
    // Check magic number. 
    tpr.read((char*)&nVal, sizeof(int));

    if (nVal != nExpectedMagicNumber)
    {
        nVal = __bswap_32(nVal);

        if (nVal == nExpectedMagicNumber)
        {
            bSwap = true;
        }
        else
        {
            std::cerr << "Could not read magic number. Aborting" << std::endl;
        }
    }
    
    tpr.read((char*)&nVal, sizeof(int));
    if (bSwap) nVal = __bswap_32(nVal);

    tpr.seekp(nInitStepLoc, std::ios_base::beg);
    nVal = 0;
    tpr.write((char*)&nVal, sizeof(int));

    tpr.seekp(nAdditionalStepLoc, std::ios_base::beg);
    tpr.read((char*)&nVal, sizeof(int));
    if (bSwap) nVal = __bswap_32(nVal);
    tpr.seekp(-sizeof(int), std::ios_base::cur);
    if (bSwap) nAdditionalSteps = __bswap_32(nAdditionalSteps);
    tpr.write((char*)&nAdditionalSteps, sizeof(int));

    tpr.seekp(nPullValueLoc, std::ios_base::beg);
    int nPullValue = *(reinterpret_cast<int*>(&fPullValue));
    if (bSwap) nPullValue = __bswap_32(nPullValue);
    tpr.write((char*)&nPullValue, sizeof(int));

    tpr.seekp(nPullKValueLoc, std::ios_base::beg);
    int nPullKValue = *(reinterpret_cast<int*>(&fPullKValue));
    if (bSwap) nPullKValue = __bswap_32(nPullKValue);
    tpr.write((char*)&nPullKValue, sizeof(int));

    tpr.close();

}


/***************************************************************************//*!
*//****************************************************************************/
int 
main(int argc, char* argv[])
{
    if (argc != 9)
    {
        std::cerr << "Error: Usage " << argv[0] 
            << " tprFile additionalStepLoc additionalSteps pullValueLoc pullValue" 
            << " pullKValueLoc pullKValue initStepLoc " << std::endl;
        exit(-1);
    }

    int nAdditionalStepLoc = atoi(argv[2]);
    int nAdditionalSteps = atoi(argv[3]);
    
    if (!nAdditionalStepLoc || nAdditionalStepLoc == INT_MAX 
        || nAdditionalStepLoc == INT_MIN)
    {
        std::cerr << "Error: Could not parse " << argv[2] 
            << " into non-zero value in range" << std::endl;
        exit(-1);
    }

    if (!nAdditionalSteps || nAdditionalSteps == INT_MAX 
        || nAdditionalSteps == INT_MIN)
    {
        std::cerr << "Error: Could not parse " << argv[3] 
            << " into non-zero value in range" << std::endl;
        exit(-1);
    }

    int nPullValueLoc = atoi(argv[4]);

    if (!nPullValueLoc || nPullValueLoc == INT_MAX 
        || nPullValueLoc == INT_MIN)
    {
        std::cerr << "Error: Could not parse " << argv[4] 
            << " into non-zero value in range" << std::endl;
        exit(-1);
    }

    float fPullValue = (float)atof(argv[5]);

    if (fPullValue == HUGE_VAL)
    {
        std::cerr << "Error: Could not parse " << argv[5] << " into a float." 
            << std::endl;
        exit(-1);
    }

    int nPullKValueLoc = atoi(argv[6]);

    if (!nPullKValueLoc || nPullKValueLoc == INT_MAX 
        || nPullKValueLoc == INT_MIN)
    {
        std::cerr << "Error: Could not parse " << argv[6] 
            << " into non-zero value in range" << std::endl;
        exit(-1);
    }

    float fPullKValue = (float)atof(argv[7]);

    if (fPullKValue == HUGE_VAL)
    {
        std::cerr << "Error: Could not parse " << argv[7] << " into a float." 
            << std::endl;
        exit(-1);
    }

    int nInitStepLoc = atoi(argv[8]);

    updateTPR(argv[1], nAdditionalStepLoc, nAdditionalSteps, nPullValueLoc, 
        fPullValue, nPullKValueLoc, fPullKValue, nInitStepLoc);
}
