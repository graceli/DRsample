#include <byteswap.h>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits.h>

static const unsigned int nExpectedMagicNumber = 14;
static const int tpx_generation = 17;

enum {
  F_BONDS,
  F_G96BONDS,
  F_MORSE,
  F_CUBICBONDS,
  F_CONNBONDS,
  F_HARMONIC,
  F_FENEBONDS,
  F_TABBONDS,
  F_TABBONDSNC,
  F_ANGLES, 
  F_G96ANGLES,
  F_CROSS_BOND_BONDS,
  F_CROSS_BOND_ANGLES,
  F_UREY_BRADLEY,
  F_QUARTIC_ANGLES,
  F_TABANGLES,
  F_PDIHS,
  F_RBDIHS, 
  F_FOURDIHS,
  F_IDIHS, 
  F_PIDIHS, 
  F_TABDIHS,
  F_LJ14,
  F_COUL14,
  F_LJC14_Q,
  F_LJC_PAIRS_NB,
  F_LJ,
  F_BHAM,
  F_LJ_LR,
  F_BHAM_LR,
  F_DISPCORR,
  F_COUL_SR,
  F_COUL_LR,
  F_RF_EXCL,
  F_COUL_RECIP,
  F_DPD,
  F_POLARIZATION,
  F_WATER_POL,
  F_THOLE_POL,
  F_POSRES,
  F_DISRES,
  F_DISRESVIOL,
  F_ORIRES,
  F_ORIRESDEV,
  F_ANGRES,
  F_ANGRESZ,
  F_DIHRES,
  F_DIHRESVIOL,
  F_CONSTR,
  F_CONSTRNC,
  F_SETTLE,
  F_VSITE2,
  F_VSITE3,
  F_VSITE3FD,
  F_VSITE3FAD,
  F_VSITE3OUT,
  F_VSITE4FD,
  F_VSITE4FDN,
  F_VSITEN,
  F_COM_PULL,
  F_EQM,
  F_EPOT,
  F_EKIN,
  F_ETOT,
  F_ECONSERVED,
  F_TEMP,
  F_PRES,
  F_DVDL,
  F_DKDL,
  F_DGDL_CON,
  F_NRE		/* This number is for the total number of energies	*/
};

enum {
  egcTC,    egcENER,   egcACC, egcFREEZE, 
  egcUser1, egcUser2,  egcVCM, egcXTC,
  egcORFIT, egcQMMM,
  egcNR 
};

/* Struct used to maintain tpx compatibility when function types are added */
typedef struct {
  int fvnr; /* file version number in which the function type first appeared */
  int ftype; /* function type */
} t_ftupd;

static const t_ftupd ftupd[] = {
  { 20, F_CUBICBONDS        },
  { 20, F_CONNBONDS         },
  { 20, F_HARMONIC          },
  { 34, F_FENEBONDS         },
  { 43, F_TABBONDS          },
  { 43, F_TABBONDSNC        },
  { 30, F_CROSS_BOND_BONDS  },
  { 30, F_CROSS_BOND_ANGLES },
  { 30, F_UREY_BRADLEY      },
  { 34, F_QUARTIC_ANGLES    },
  { 43, F_TABANGLES         },
  { 26, F_FOURDIHS          },
  { 26, F_PIDIHS            },
  { 43, F_TABDIHS           },
  { 41, F_LJC14_Q           },
  { 41, F_LJC_PAIRS_NB      },
  { 32, F_BHAM_LR           },
  { 32, F_RF_EXCL           },
  { 32, F_COUL_RECIP        },
  { 46, F_DPD               },
  { 30, F_POLARIZATION      },
  { 36, F_THOLE_POL         },
  { 22, F_DISRESVIOL        },
  { 22, F_ORIRES            },
  { 22, F_ORIRESDEV         },
  { 26, F_DIHRES            },
  { 26, F_DIHRESVIOL        },
  { 49, F_VSITE4FDN         },
  { 50, F_VSITEN            },
  { 46, F_COM_PULL          },
  { 20, F_EQM               },
  { 46, F_ECONSERVED        },
  { 54, F_DGDL_CON          }
};

#define asize(a) (sizeof(a)/sizeof((a)[0]))
#define NFTUPD asize(ftupd)
#define MAXNODES 256

/***************************************************************************//*!
*//****************************************************************************/
void
doFFParams(std::fstream& tpr, int nFileVersion, int nRealSize, bool bSwap)
{
    int nPos = 0;
    int nVal = 0;
    tpr.seekp(sizeof(int), std::ios_base::cur);
    if (nFileVersion < 57)
    {
        tpr.seekp(sizeof(int), std::ios_base::cur);
    }
    int nFFTypes;
    tpr.read((char*)&nFFTypes, sizeof(int)); 
    if (bSwap) nFFTypes = __bswap_32(nFFTypes);
    nPos = 0;
    for (int i=0; i != nFFTypes; ++i)
    {
        tpr.read((char*)&nVal, sizeof(int)); 
        if (bSwap) nVal = __bswap_32(nVal);
        switch(nVal)
        {
            case F_ANGLES:
            case F_G96ANGLES:
            case F_BONDS:
            case F_G96BONDS:
            case F_HARMONIC:
            case F_IDIHS:
                nPos += nRealSize * 4;
                break;
            case F_FENEBONDS:
                nPos += nRealSize * 2;
                break;
            case F_TABBONDS:
            case F_TABBONDSNC:
            case F_TABANGLES:
            case F_TABDIHS:
                nPos += nRealSize * 2 + sizeof(int);
                break;
            case F_CROSS_BOND_BONDS:
                nPos += nRealSize * 3;
                break;
            case F_CROSS_BOND_ANGLES:
            case F_UREY_BRADLEY:
                nPos += nRealSize * 4;
                break;
            case F_QUARTIC_ANGLES:
                nPos += nRealSize * 6; ///////
                break;
            case F_BHAM:
            case F_MORSE:
            case F_CUBICBONDS:
                nPos += nRealSize * 3;
                break;
            case F_CONNBONDS:
                break;
            case F_POLARIZATION:
                nPos += nRealSize;
                break;
            case F_WATER_POL:
            if (nFileVersion < 31) 
            {
                std::cerr << "Unsupported TPR. Check tpxio.c" << std::endl;
                exit(-1);
            }
                nPos += nRealSize * 6;
                break;
            case F_THOLE_POL:
                nPos += nRealSize * 4;
                break;
            case F_LJ:
                nPos += nRealSize * 2;
                break;
            case F_LJ14:
                nPos += nRealSize * 4;
                break;
            case F_LJC14_Q:
                nPos += nRealSize * 5;
                break;
            case F_LJC_PAIRS_NB:
                nPos += nRealSize * 4;
                break;
            case F_PDIHS:
            case F_PIDIHS:
            case F_ANGRES:
            case F_ANGRESZ:
                nPos += nRealSize * 2;
                if ((nVal == F_ANGRES || nVal == F_ANGRESZ) && nFileVersion < 42) 
                {
                    nPos += nRealSize * 2;
                }
                else 
                {
                    nPos += nRealSize * 2 + sizeof(int);
                }
                break;
            case F_DISRES:
                nPos += nRealSize * 4 + sizeof(int) * 2;
                break;
            case F_ORIRES:
                nPos += nRealSize * 3 + sizeof(int) * 3;
                break;
            case F_DIHRES:
                nPos += nRealSize * 3 + sizeof(int) * 2;
                break;
            case F_POSRES:
                nPos += nRealSize * 6;
                if (nFileVersion >= 27) 
                {
                    nPos += nRealSize * 6;
                }
                break;
            case F_RBDIHS:
                nPos += nRealSize * 6;
                if(nFileVersion >= 25) 
                    nPos += nRealSize * 6;
                break;
            case F_FOURDIHS:
                nPos += nRealSize * 12;
                break;
            case F_CONSTR:
            case F_CONSTRNC:
            case F_SETTLE:
                nPos += nRealSize * 2;
                break;
            case F_VSITE2:
                nPos += nRealSize;
                break;
            case F_VSITE3:
            case F_VSITE3FD:
            case F_VSITE3FAD:
                nPos += nRealSize * 2;
                break;
            case F_VSITE3OUT:
            case F_VSITE4FD: 
            case F_VSITE4FDN: 
                nPos += nRealSize * 3;
                break;
            case F_VSITEN:
                nPos += nRealSize + sizeof(int);
                break;
            default:
                std::cout << "Unhandled type: " << nVal << std::endl;
        }
    }
    if (nFileVersion >= 57)
    {
        tpr.seekp(nRealSize, std::ios_base::cur);
    }
    
    tpr.seekp(nPos, std::ios_base::cur);
}


/***************************************************************************//*!
*//****************************************************************************/
void
doBlock(std::fstream& tpr, int nFileVersion, bool bSwap)
{
    if (nFileVersion < 44)
    {
        tpr.seekp(sizeof(int)*MAXNODES, std::ios_base::cur);
    }

    int blocknr;
    tpr.read((char*)&blocknr, sizeof(int));
    if (bSwap) blocknr = __bswap_32(blocknr);

    int dummy;
    if (nFileVersion < 51)
    {
        tpr.read((char*)&dummy, sizeof(int));
        if (bSwap) dummy = __bswap_32(dummy);
    }

    tpr.seekp(sizeof(int) * (blocknr+1), std::ios_base::cur);

    if (nFileVersion < 51)
    {
        tpr.seekp(sizeof(int) * dummy, std::ios_base::cur);
    }
}

/***************************************************************************//*!
*//****************************************************************************/
void
doBlockA(std::fstream& tpr, int nFileVersion, bool bSwap)
{
    if (nFileVersion < 44)
    {
        tpr.seekp(sizeof(int) * MAXNODES, std::ios_base::cur);
    }

    int nBlockNR, nBlockNRA;
        
    tpr.read((char*)&nBlockNR, sizeof(int)); 
    if (bSwap) nBlockNR = __bswap_32(nBlockNR);

    tpr.read((char*)&nBlockNRA, sizeof(int)); 
    if (bSwap) nBlockNRA = __bswap_32(nBlockNRA);

    tpr.seekp(sizeof(int) * (nBlockNR+1), std::ios_base::cur);
    tpr.seekp(sizeof(int) * nBlockNRA, std::ios_base::cur);
}



/***************************************************************************//*!
*//****************************************************************************/
void
doInputRec(std::fstream& tpr, int nFileVersion, int nRealSize, bool bSwap)
{
    // Assuming that fileversion >= 1.
    int nVal;
    tpr.seekp(sizeof(int), std::ios_base::cur);
    int nSteps;
    int nStepPos;
    int nPullPos;
    int nPullK1Pos;
    int nInitStepPos = 0;
    float fPull;
    nStepPos = tpr.tellp();
    tpr.read((char*)&nSteps, sizeof(int)); 
    if (bSwap) nSteps = __bswap_32(nSteps); // OK to here. 

    int nSeek = 0;
    if (nFileVersion > 25)
    {
        // init_step
        nInitStepPos = tpr.tellp();
        nSeek += sizeof(int);
    }

    if (nFileVersion >= 58)
    {
        // simulation_part
        nSeek += sizeof(int);
    }

    if (nFileVersion < 53)
    {
        nSeek += sizeof(int);
        if (nFileVersion >= 45)
        {
            nSeek += sizeof(int);
        }
    }

    // ns_type, nstlist, ndelta
    nSeek += sizeof(int) * 3;

    if (nFileVersion < 41)
    {
        nSeek += sizeof(int) * 2;
    }

    if (nFileVersion >= 45)
    {
        nSeek += nRealSize;
    }

    // nstcomm
    nSeek += sizeof(int);
    if (nFileVersion > 34)
    {
        // comm_mode
        nSeek += sizeof(int);
    }

    if (nFileVersion > 25)
    {
        // checkpoint
        nSeek += sizeof(int);
    }

    // nstcgsteep
    nSeek += sizeof(int);
    if (nFileVersion > 30)
    {
        //nbfgscorr
        nSeek += sizeof(int);
    }

    nSeek += sizeof(int) * 6 + nRealSize * 3;

    if (nFileVersion < 19)
    {
        nSeek += sizeof(int) * 2;
    }
      
    if (nFileVersion < 18)
    {
        //idum
        nSeek += sizeof(int);
    }

    // rlist, coulombtype
    nSeek += sizeof(int) + nRealSize;

    nSeek += sizeof(int) * 2 + nRealSize * 5;

    if (nFileVersion >= 37)
    {
        // epsilon_rf
        nSeek += nRealSize;
    }
    
    if (nFileVersion >= 29)
    {
        nSeek += nRealSize;
    }

    if (nFileVersion > 25)
    {
        nSeek += sizeof(int) * 3 + nRealSize * 2;
    }

    if (nFileVersion >= 55)
    {
        nSeek += nRealSize * 5;
    }

    nSeek += sizeof(int) * 4 + nRealSize;

    if (nFileVersion >= 24)
    {
        // ewald_geometry
        nSeek += sizeof(int);
    }

    if (nFileVersion <= 17)
    {
        if (nFileVersion == 17)
        {
            nSeek += sizeof(int);
        }
    }
    else
    {
        // Epsilon surface
        nSeek += nRealSize;
    }

    // boptfft ... etc
    nSeek += sizeof(int) * 3; 
    

    if (nFileVersion <= 15)
    {
        nSeek += sizeof(int);
    }

    if (nFileVersion <= 17)
    {
        nSeek += sizeof(int);
        if (nFileVersion <= 15)
        {
            nSeek += sizeof(int);
        }
    }
    else
    {
        // epc, epct
        nSeek += sizeof(int) * 2;
    }

    // Checkpoint.
    tpr.seekp(nSeek, std::ios_base::cur);
    nSeek = 0;
    
    tpr.read((char*)&nVal, sizeof(int));
    if (bSwap) nVal = __bswap_32(nVal);

    // Looks ok to here. 
    if (nFileVersion <= 15)
    {
        nSeek += nRealSize * 3 * 2;
    }
    else
    {
        nSeek += nRealSize * 3 * 3 * 2; // ref_p and compress
    }

    if (nFileVersion >= 47)
    {
        nSeek += sizeof(int) + nRealSize * 3 * 2;
    }

    if (nFileVersion > 25)
    {
        nSeek += sizeof(int); // andersen_seed
    }

    if (nFileVersion < 26)
    {
        nSeek += sizeof(int) + nRealSize; // bSimAnn, zerotemptime.
    }

    if (nFileVersion < 37)
    {
        nSeek += nRealSize; //rdum
    }

    nSeek += nRealSize; //shake_tol
    
    if (nFileVersion < 54)
    {
        nSeek += nRealSize; //fudgeQQ
    }

    nSeek += sizeof(int) + nRealSize * 2; // efep. {init,delta}_lambda

    if (nFileVersion >= 13)
    {
        nSeek += nRealSize;
    }

    if (nFileVersion >= 38)
    {
        nSeek += sizeof(int);
    }

    if (nFileVersion >= 15)
    {
        nSeek += nRealSize;
    }

    if (nFileVersion >= 57)
    {
        nSeek += sizeof(int); // eDisre
    }

    nSeek += sizeof(int); // eDisreWeighting

    nSeek += sizeof(int) * 2 + nRealSize * 2;// bDisreMixed .. nstdisreout

    if (nFileVersion >= 22)
    {
        nSeek += sizeof(int) + nRealSize * 2;
    }

    if (nFileVersion >= 26)
    {
        nSeek += nRealSize;
        if (nFileVersion < 56)
        {
             nSeek += sizeof(int) + nRealSize;
        }
    }

    nSeek += nRealSize * 2; // em_stepsize em_tol

    if (nFileVersion >= 22)
    {
         nSeek += sizeof(int);
    }
    if (nFileVersion >= 11)
    {
         nSeek += sizeof(int); //niter
    }
      
    if (nFileVersion >= 21)
    {
         nSeek += nRealSize;
    }

    nSeek += sizeof(int) * 2 + nRealSize; //eConstrAlg .. LincsWarnAngle

    if (nFileVersion <= 14)
    {
        nSeek += sizeof(int);
    }

    if (nFileVersion >= 26)
    {
        nSeek += sizeof(int);
    }

    if (nFileVersion < 33)
    {
        nSeek += nRealSize; // bd_temp
    }
    nSeek += sizeof(int) + nRealSize; //bd_fric ld_seed

    if (nFileVersion >= 33)
    {
        nSeek += nRealSize * 3 * 3; //deform
    }

    if (nFileVersion >= 14)
    {
        nSeek += nRealSize; // cos_accel;
    }
        
    nSeek += sizeof(int) * 4 + nRealSize * 4;

    tpr.seekp(nSeek, std::ios_base::cur);
    nSeek = 0;

    int ePull;
    tpr.read((char*)&ePull, sizeof(int));
    if (bSwap) ePull = __bswap_32(ePull);

    // Pull code.
    int ngrp;
    tpr.read((char*)&ngrp, sizeof(int));
    if (bSwap) ngrp = __bswap_32(ngrp);
    tpr.seekp(sizeof(int) * 6 + nRealSize * 3, std::ios_base::cur);

    for (int g=0; g != ngrp+1; ++g)
    {
        int nat;
        tpr.read((char*)&nat, sizeof(int));
        if (bSwap) nat = __bswap_32(nat);

        tpr.seekp(sizeof(int) * nat, std::ios_base::cur);

        int weight;
        tpr.read((char*)&weight, sizeof(int));
        if (bSwap) weight = __bswap_32(weight);

        tpr.seekp(nRealSize * weight, std::ios_base::cur);

        // pbcatom, vec
        tpr.seekp(sizeof(int) + nRealSize * 3, std::ios_base::cur);

        tpr.seekp(nRealSize * 3, std::ios_base::cur);

        if (g == 1)
        {
            nPullPos = tpr.tellp();
            nPullPos -= nRealSize;
            tpr.seekp(-nRealSize, std::ios_base::cur);

            if (nRealSize == 4)
            {
                int nPull;
                tpr.read((char*)&nPull, sizeof(int));
                if (bSwap) nPull = __bswap_32(nPull);

                float fPull = *(reinterpret_cast<float*>(&nPull));
//                std::cout << "Pull " << fPull << std::endl;
            }
            else
            {
                std::cerr << "Last of the double precision not implemented yet." << std::endl;
                exit(-1);
            }

            tpr.seekp(nRealSize, std::ios_base::cur);

            nPullK1Pos = tpr.tellp();

            if (nRealSize == 4)
            {
                int nPullK1;
                tpr.read((char*)&nPullK1, sizeof(int));
                if (bSwap) nPullK1 = __bswap_32(nPullK1);

                float fPullK1 = *(reinterpret_cast<float*>(&nPullK1));
//                std::cout << "Pull_K1 " << fPullK1 << std::endl;
            }
            else
            {
                std::cerr << "Last of the double precision not implemented yet." << std::endl;
                exit(-1);
            }
        }
        else
        {
            // Skip over rate and k.
            tpr.seekp(nRealSize * 2, std::ios_base::cur);
        }

        if (nFileVersion >= 56)
        {
            tpr.seekp(nRealSize, std::ios_base::cur);
        }
    }

    std::cout << "Step " << nStepPos << std::endl;
    std::cout << "Pull " << nPullPos << std::endl;
    std::cout << "Pull_K1 " << nPullK1Pos << std::endl;
    std::cout << "InitStep " << nInitStepPos << std::endl;


}


/***************************************************************************//*!
*//****************************************************************************/
void
updateTPR(char* sFile)
{
    char buf[1024];
    bool bSwap = false;
    int nVal, nRealSize, nFileVersion, nFileGeneration;
    int nBox, ngtc, ntop, nNumAtoms, nInputRec, nX, nV, nF, nStep;
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

    tpr.read(buf, nVal); // Read gromacs version.
    int nPos = tpr.tellp();
    if (nPos % 4)
    {
        nPos += 4 - (nPos % 4); // 32-bit align
        tpr.seekp(nPos);
    }

    // Looks like some tpr files may have 64-bit reals, and others
    // will have 32-bit ones.
    tpr.read((char*)&nRealSize, sizeof(int)); 
    if (bSwap) nRealSize = __bswap_32(nRealSize);

    tpr.read((char*)&nFileVersion, sizeof(int));
    if (bSwap) nFileVersion = __bswap_32(nFileVersion);

    if (nFileVersion >= 26)
    {
        tpr.read((char*)&nFileGeneration, sizeof(int));
        if (bSwap) nFileGeneration = __bswap_32(nFileGeneration);
    }
    else
    {
        nFileGeneration = 0;
    }

    tpr.read((char*)&nNumAtoms, sizeof(int)); 
    if (bSwap) nNumAtoms = __bswap_32(nNumAtoms);


    if (nFileVersion >= 28)
    {
        tpr.read((char*)&ngtc, sizeof(int));
        if (bSwap) ngtc = __bswap_32(ngtc);            
    }
    else
    {
        ngtc = 0;
    }

    tpr.read((char*)&nStep, sizeof(int)); 
    if (bSwap) nStep = __bswap_32(nStep);

    tpr.seekp(nRealSize * 2, std::ios_base::cur);

    tpr.read((char*)&nInputRec, sizeof(int)); 
    if (bSwap) nInputRec = __bswap_32(nInputRec);    

    tpr.read((char*)&ntop, sizeof(int)); 
    if (bSwap) ntop = __bswap_32(ntop);

    tpr.read((char*)&nX, sizeof(int)); 
    if (bSwap) nX = __bswap_32(nX);

    tpr.read((char*)&nV, sizeof(int)); 
    if (bSwap) nV = __bswap_32(nV);

    tpr.read((char*)&nF, sizeof(int)); 
    if (bSwap) nF = __bswap_32(nF);    

    tpr.read((char*)&nBox, sizeof(int)); 
    if (bSwap) nBox = __bswap_32(nBox);

    // fp is at the end of the header right now. 
    nPos = 0;
    if (nBox)
    {
        nPos += nRealSize * 9; // Box vectors.
        if (nFileVersion >= 51)
        {
            nPos += nRealSize * 9; // box_rel.
        }
        if (nFileVersion >= 28)
        {
            nPos += nRealSize * 9; // boxv
            if (nFileVersion < 56)
            {
                nPos += nRealSize * 9;
            }
        }
    }

    if (ngtc > 0 && nFileVersion >= 28)
    {
        nPos += ngtc * nRealSize * 2; // nosehoover & tcouple.
    }

    if (nFileVersion < 26)
    {
        std::cerr << "Error. File version < 26. This hasn't been implemented yet. " 
            << std::endl;
        exit(-1);
    }

    tpr.seekp(nPos, std::ios_base::cur);

    // Read symbol table.
    if (ntop)
    {
        tpr.read((char*)&nVal, sizeof(int)); 
        if (bSwap) nVal = __bswap_32(nVal);

        int nSymtabSize = nVal;

        for (int i=0; i != nSymtabSize; ++i)
        {
            tpr.seekp(sizeof(int), std::ios_base::cur);
            tpr.read((char*)&nVal, sizeof(int)); 
            if (bSwap) nVal = __bswap_32(nVal);
            if (nVal % 4)
            {
                nVal += 4 - (nVal % 4);
            }
            tpr.seekp(nVal, std::ios_base::cur);
        }
    }

    
    tpr.seekp(sizeof(int), std::ios_base::cur);

    if (nFileVersion >= 57)
    {
        doFFParams(tpr, nFileVersion, nRealSize, bSwap);
    }
    
//    tpr.seekp(nPos, std::ios_base::cur);

    int nmoltype;
    tpr.read((char*)&nmoltype, sizeof(int)); 
    if (bSwap) nmoltype = __bswap_32(nmoltype); // fp now done reading moltype.

    for (int i=0; i != nmoltype; ++i) // do_moltype from do_mtop
    {
        int nr, nres, ngrpname;
        if (nFileVersion >= 57)
        {
            tpr.seekp(sizeof(int), std::ios_base::cur); 
        }

        // Start of do_atoms from do_moltype
        tpr.read((char*)&nr, sizeof(int)); 
        if (bSwap) nr = __bswap_32(nr);
        
        tpr.read((char*)&nres, sizeof(int)); 
        if (bSwap) nres = __bswap_32(nres); // OK to here. 

        if (nFileVersion < 57)
        {
            tpr.read((char*)&ngrpname, sizeof(int)); 
            if (bSwap) ngrpname = __bswap_32(ngrpname);
        }

        nPos = tpr.tellp();

        // This can be condensed iff the atom info isn't required later on.
        for (int j=0; j != nr; ++j)
        {
            // Note - do_uchar etc still takes 4 bytes on the disk!.
            int nCount = nRealSize * 4 + sizeof(int) * 1 + sizeof(unsigned int) * 2;
            tpr.seekp(nCount, std::ios_base::cur); 
        
            int resnr;
            tpr.read((char*)&resnr, sizeof(int)); 
            if (bSwap) resnr = __bswap_32(resnr); // OK to here. 

            int myngrp = egcNR;
            if (nFileVersion < 23) 
                myngrp = 8;
            else if (nFileVersion < 39) 
                myngrp = 9;
    
            if (nFileVersion >= 52)
            {
                int atomnr;
                tpr.read((char*)&atomnr, sizeof(int)); 
                if (bSwap) atomnr = __bswap_32(atomnr); 
            }

            nCount = 0;
            if (nFileVersion < 57)
            {
                nCount += sizeof(unsigned char) * myngrp;
            }
            tpr.seekp(nCount, std::ios_base::cur); // end of do_atom.
        }

        tpr.seekp(sizeof(int) * nr, std::ios_base::cur);

        if (nFileVersion > 20)
        {
            tpr.seekp(sizeof(int) * 2 * nr, std::ios_base::cur); // 2 do_strstr
        }

        tpr.seekp(sizeof(int) * nres, std::ios_base::cur); //do_strstr

        if (nFileVersion < 57)
        {
            tpr.seekp(sizeof(int) * ngrpname, std::ios_base::cur); // do_strstr end of do_atoms

            int myngrp; // Starting do_grps from do_atoms
            int ngrp = egcNR;
            if (nFileVersion < 23) 
                myngrp = 8;
            else if (nFileVersion < 39) 
                myngrp = 9;
            else
                myngrp = ngrp;
  
            for(int j=0; (j<ngrp); j++) 
            {
                if (j<myngrp) 
                {
                    tpr.read((char*)&nVal, sizeof(int)); 
                    if (bSwap) nVal = __bswap_32(nVal);
                    tpr.seekp(nVal * sizeof(int), std::ios_base::cur);
                }
            }
        } // Done do_atoms in do_moltype

        if (nFileVersion >= 57)
        {
            for(int j=0; (j<F_NRE); j++)  // do_ilists
            {
                int k;
                bool bClear = false;
                for (k=0; k<NFTUPD; k++)
                {
                	if ((nFileVersion < ftupd[k].fvnr) && (j == ftupd[k].ftype))
                    {
                        bClear = true;
                    }
                }

                if (!bClear)
                {
                    tpr.read((char*)&nVal, sizeof(int)); 
                    if (bSwap) nVal = __bswap_32(nVal);
                    tpr.seekp(sizeof(int) * nVal, std::ios_base::cur);
                }
            } // End do_ilists

            doBlock(tpr, nFileVersion, bSwap);
        }

        // do_blocka
        doBlockA(tpr, nFileVersion, bSwap);
    } // OK to here. End of do_moltype loop.

    int nmolblock;
    if (nFileVersion >= 57)
    {
        tpr.read((char*)&nmolblock, sizeof(int)); 
        if (bSwap) nmolblock = __bswap_32(nmolblock);
    }
    else
    {
        nmolblock = 1;
    }

    if (nFileVersion >= 57)
    {
        for (int i=0; i != nmolblock; ++i)
        {
            // do_molblock
            tpr.seekp(sizeof(int)*3, std::ios_base::cur); //type, nmol, natoms_mol
            int xA;
    
            tpr.read((char*)&xA, sizeof(int)); 
            if (bSwap) xA = __bswap_32(xA);

            if (xA > 0)
            {
                tpr.seekp(nRealSize * 3 * xA, std::ios_base::cur);
            }
            
            int xB;
    
            tpr.read((char*)&xB, sizeof(int)); 
            if (bSwap) xB = __bswap_32(xB);

            if (xB > 0)
            {
                tpr.seekp(nRealSize * 3 * xB, std::ios_base::cur);
            }
        } // End do_molblock
        tpr.seekp(sizeof(int), std::ios_base::cur); // mtop->natoms
    }

    // do_atomtypes
    if (nFileVersion > 25)
    {
        tpr.read((char*)&nVal, sizeof(int)); 
        if (bSwap) nVal = __bswap_32(nVal);

        tpr.seekp(nVal * 3 * nRealSize, std::ios_base::cur);

        if (nFileVersion >= 40)
        {
            tpr.seekp(sizeof(int) * nVal, std::ios_base::cur);
        }
    }

    if (nFileVersion < 57) // do_idef
    {
        doFFParams(tpr, nFileVersion, nRealSize, bSwap);
        if (nFileVersion >= 54)
        {
            tpr.seekp(nRealSize, std::ios_base::cur);
        }
        
        for(int j=0; (j<F_NRE); j++)  // do_ilists
        {
            int k;
            bool bClear = false;
            for (k=0; k<NFTUPD; k++)
            {
            	if ((nFileVersion < ftupd[k].fvnr) && (j == ftupd[k].ftype))
                {
                    bClear = true;
                }
            }

            if (!bClear)
            {
                tpr.read((char*)&nVal, sizeof(int)); 
                if (bSwap) nVal = __bswap_32(nVal);
                tpr.seekp(sizeof(int) * nVal, std::ios_base::cur);
            }
        } // End do_ilists        
    } // End do_idef

    if (nFileVersion >= 57)
    {
        // do_groups
        int myngrp; // Starting do_grps from do_atoms
        int ngrp = egcNR;
        if (nFileVersion < 23) 
            myngrp = 8;
        else if (nFileVersion < 39) 
            myngrp = 9;
        else
            myngrp = ngrp;
  
        for(int j=0; (j<ngrp); j++) 
        {
            if (j<myngrp) 
            {
                tpr.read((char*)&nVal, sizeof(int)); 
                if (bSwap) nVal = __bswap_32(nVal);
                tpr.seekp(nVal * sizeof(int), std::ios_base::cur);
            }
        } // end do_grps

        int ngrpname;
        tpr.read((char*)&ngrpname, sizeof(int)); 
        if (bSwap) ngrpname = __bswap_32(ngrpname); // OK to here. 

        tpr.seekp(sizeof(int) * ngrpname, std::ios_base::cur);
        for (int g=0; g != egcNR; ++g)
        {
            int ngrpnr;
            tpr.read((char*)&ngrpnr, sizeof(int)); 
            if (bSwap) ngrpnr = __bswap_32(ngrpnr); 

            if (ngrpnr)
            {
                // Again, this is a char, but look like it gets 
                // promoted to int before storing to disk. 
                // Should confirm this sometime.
                tpr.seekp(sizeof(int) * ngrpnr, std::ios_base::cur); 
            }
        } // OK to here...?
        
    } // End do_groups

    if (nFileVersion < 57)
    {
        doBlock(tpr, nFileVersion, bSwap); 
        doBlock(tpr, nFileVersion, bSwap); 
    }

    if (nFileVersion < 51)
    {
        doBlockA(tpr, nFileVersion, bSwap);
    }

    // Done do_mtop now.
    if (nX)
    {
        tpr.seekp(nRealSize * nNumAtoms * 3, std::ios_base::cur);
    }

    if (nV)
    {
        tpr.seekp(nRealSize * nNumAtoms * 3, std::ios_base::cur);
    }

    if (nF)
    {
        tpr.seekp(nRealSize * nNumAtoms * 3, std::ios_base::cur);
    }


    if (nFileVersion >= 26)
    {
        if (nInputRec)
        {
            if (nFileVersion >= 53)
            {
                int epbc, period;
                tpr.read((char*)&epbc, sizeof(int)); 
                if (bSwap) epbc = __bswap_32(epbc); 

                tpr.read((char*)&period, sizeof(int)); 
                if (bSwap) period = __bswap_32(period); 
                
            }

            if (nFileGeneration <=  tpx_generation)
            {
                doInputRec(tpr, nFileVersion, nRealSize, bSwap);
            }
        }
    }

    tpr.close();
}


/***************************************************************************//*!
*//****************************************************************************/
int 
main(int argc, char *argv[])
{
    if (argc != 2)
    {
        std::cerr << "Error: Usage " << argv[0] 
            << " tprFile " << std::endl;
        exit(-1);
    }

    updateTPR(argv[1]);

}
