/*
 *  tempest.h
 *
 *  Created by Brendan Faherty on 12/8/08.
 *  Copyright 2008 Dartmouth College. All rights reserved.
 *
 */

#ifndef TEMPEST_H
#define TEMPEST_H

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdio.h>
#include <cstdlib>
#include <fcntl.h>
#include <ctype.h>
#include <errno.h>
#include <string.h>
#include <cmath>
#include <float.h>
#include <time.h>
#include <vector>
#include <stack>
#include <map>

typedef off_t f_off; //MSToolkit 
#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"

#ifdef _WIN32
#include <io.h>
//typedef size_t off_t;
typedef char bool;
#define isnan(n) _isnan(n)
#define strcasecmp _stricmp
#define read _read
#define NULL_FILE_HANDLE "nul"
#else
#include <strings.h>
#include <stdbool.h>
#define NULL_FILE_HANDLE "/dev/null"
#endif

#define VERSION_STRING "1.0"
#define LICENSE_STRING "(c) 2009 Dartmouth College"
#define AUTHORS_STRING "Jeffrey Milloy, Brendan Faherty"

//==================================================================================================
//  Settings and Defaults
//==================================================================================================

#define PROFILE 0

#define STRING_SIZE 256
#define CASE_SHIFT 32

// params
#define DEFAULT_PARAMS_FILE "tempest.params"
// kernels
#define KERNEL_FILE "kernels.cl"

// chosen platform/device
#define PLATFORM_ID 0
#define DEVICE_ID 0

// msms
#define MAX_CHARGE 4
#define INITIAL_PEAKS_BUFFER 100
#define MIN_DIGEST_MASS 400 //TODO remove/rename
#define MAX_DIGEST_MASS 5000 //TODO remove/rename
#define MAX_DIGEST_MASS_MS2 6000
#define MSMS_NORM_REGIONS 50
#define MSMS_NOISE_CUTOFF 0.0
#define CROSS_CORRELATION_WINDOW 75

// digestion
#define MAX_LENGTH_REFERENCE 64
#define MIN_PEPTIDE_LENGTH 1
#define MAX_PEPTIDE_LENGTH 64
#define MAX_MOD_TARGETS_PER_PEPTIDE 15
#define INITIAL_PROTEIN_BUFFER 2048

// scoring
#define PRIMARY_INTENSITY 50
#define PHOS_NL_INTENSITY 50
#define FLANKING_INTENSITY 25
#define NL_INTENSITY 10
#define MIN_SCORE -1000.0f

// config
extern size_t MAX_BLOCKDIM;
extern size_t DEFAULT_CANDIDATE_BUFFER_SIZE;
extern size_t DEFAULT_BLOCKDIM_SCORE;
extern size_t BLOCKDIM_PREBUILD;
extern size_t BLOCKDIM_PREBUILD_TRANSFORMATION;
extern size_t BLOCKDIM_BUILD;
extern size_t BLOCKDIM_TRANSFORM;
extern size_t BLOCKDIM_REDUCE;
extern size_t CL_ONE;

// constants
#define PROTON_MASS     1.00727642
#define H_MASS          1.007825032
#define C_MASS          12.0
#define N_MASS          14.003074007
#define O_MASS          15.994914622

#define OH_MASS         17.002739654
#define NH3_MASS        17.026549103
#define H2O_MASS        18.010564686
#define CO_MASS         27.994914622
#define H3PO4_MASS      97.97689509

#define KB 1024.0f
#define MB 1048576.0f //1024*1024
#define INT_MAX 2147483647;

//==================================================================================================
// Macros
//==================================================================================================

#define aa2i(c) (int) (c - 'A')
#define safe_free(p) if (p) { free(p); p=0; }
#define safe_cudaFree(p) if (p) { cudaFree(p); p=0; }
#define safe_cudaFreeHost(p) if (p) { cudaFreeHost(p); p=0; }

//==================================================================================================
//  Structures
//==================================================================================================

enum SETUP_MODE { NO_PREBUILD, SEQ_PREBUILD, PAR_PREBUILD };
enum SCORE_MODE { SCORE, BUILD_SCORE, BUILD_TRANS_SCORE,
                  SCORE_SHARED, BUILD_SCORE_SHARED, BUILD_TRANS_SCORE_SHARED };

typedef struct MOD {
    char    cAminoAcid;
    char    cSymbol;
    double  dMassDiff;
} mod_t;

struct ARGS {
    char *sFasta;
    char *sSpectra;
    char *sOut;
    char *sParams;

    int iPrintLevel;
    bool bPrintCandidates;
    bool bSplitOutput;
    bool bNoDigest;
    bool bNoMatch;
    bool bNoScore;
};

struct PARAMS {
    // digestion
    char  *sDigestSites;
    char  *sDigestNoSites;
    int    iDigestOffset;
    int    iDigestSpecificity;
    int    iDigestMaxMissedCleavages;
    int    iMinPeptideLength;
    int    iMaxPeptideLength;

    // modifications
    mod_t *tMods;
    int    iNumMods;
    int    iModificationsMax;

    //neutral loss counts
    int   numAANL;
    int   numNtermNL;
    int   numCtermNL;

    double dPeptideNtermMass;
    double dPeptideCtermMass;
    double dProteinNtermMass;
    double dProteinCtermMass;

    double dVariableNtermMassDiff;
    double dVariableCtermMassDiff;
    char   cVariableNtermSymbol;
    char   cVariableCtermSymbol;

    // precursor (MS1) tolerance
    float  fPrecursorTolerance;
    bool   bPrecursorTolerancePPM;

    // fragment (MS2) tolerance
    float  fFragmentTolerance;
    bool   bFragmentTolerancePPM;

    // spectra processing
    bool   bRemovePrecursors;
    float *fRemoveMzRanges;
    int    iNumRemoveMzRanges;

    // similarity scoring
    bool   bCrossCorrelation;
    bool   bTrackDuplicates;
    bool   bTheoreticalFlanking;

    // output
    int    iNumOutputPSMs;
    bool   bFixDeltaScore;
};

struct CONFIG {
    int    iPlatform;
    int    iDevice;
    bool   bForceNoGPU;
    bool   bForceShared;
    bool   bForceNoPrebuild;
    size_t iCandidateBufferSize;
    size_t iScoreBlockDim;
    size_t iScoreNumBlocks;
    // scan range
    int minScan;
    int maxScan;

    enum SETUP_MODE eSetup;
    enum SCORE_MODE eScore;
};

struct INFO {
    // MS/MS scan info
    int    iNumSpectra;
    float  fPrecursorBinWidth;
    int    iNumPrecursorBins;
    float  fMinPrecursorMass;
    float  fMaxPrecursorMass;

    // MS/MS peak info
    long   lNumMS2Peaks;
    float  fMS2BinWidth;
    int    iNumMS2Bins;
    int    iCrossCorrelationWidth;

    // Digest info
    double dMinPeptideMass;
    double dMaxPeptideMass;
    int    iNumProteins;
    long   lNumPeptides;
    long   lNumPSMs;

    // Runtime info
    time_t tStart;
};

struct GPUINFO {
    // profiling requests
    bool   bTimeGPU;
    bool   bTimeKernel;
    bool   bTimeSend;
    bool   bTimeIdle;
    bool   bTimePrebuild;
    bool   bCountSyncBlocks; //counts the times when the cpu was blocked for synchronization

    // profiling results
    int    iNumScoringKernels;
    float  fTotalTimeGPU;
    float  fTotalTimeKernel;
    float  fTotalTimeSend;
    float  fTotalTimeIdle;
    float  fTotalTimePrebuild;
    int    iNumSyncBlocks;
};

struct CLINFO {
    cl_uint  iMaxComputeUnits;
    cl_uint  iAddressBits;
    cl_bool  bUnifiedMemory;
    cl_ulong lMaxWorkGroupSize;
    cl_ulong lGlobalMemSize;
    cl_ulong lMaxMemAllocSize;
};

// Candidate Peptide
typedef struct cobject {
    int    iPeptideLength;
    char   sPeptide[MAX_PEPTIDE_LENGTH];
    bool   bNtermMod;
    bool   bCtermMod;
    float  fPeptideMass;

    int    iProtein;
    char   cBefore;
    char   cAfter;
} cObj;

// Peptide-Spectrum Match
typedef struct mobject {
    int    iPeptideLength;
    char   sPeptide[MAX_PEPTIDE_LENGTH];
    bool   bNtermMod;
    bool   bCtermMod;
    float  fPeptideMass;
    float  fScore;

    int    iProtein;
    char   cBefore;
    char   cAfter;
    int    iNumOccurrences;
} mObj;

// Observed Scan
typedef struct eobject {
    long     lIndex;
    char    *sName;
    double   dPrecursorMass;
    int      iPrecursorCharge;
    cObj    *pCandidateBufferFill;
    int      iNumBufferedCandidates;
    int      iNumCandidates;
    cl_event clEventSent;
    //cudaEvent_t cudaEventSent;
    
    struct eobject* next;
} eObj;

//Mass delta and score weighting for a neutral loss event
typedef struct nlValue {
    int   numSites;
    char  aaString[64];
    float massDelta;
    float weighting;
} nlValue;

//==================================================================================================
//  External global variables
//==================================================================================================

extern struct ARGS args;
extern struct PARAMS params;
extern struct CONFIG config;
extern struct INFO tempest;
extern struct GPUINFO gpu_info;
extern struct CLINFO cl_info;
    
extern eObj **eScanIndex;
extern bool **bScanIndex;
extern char *sProteinReferences;
    
extern double dMassAA['z'-'A'+1];
extern float fMassProton[MAX_CHARGE];
extern float fMassNH3[MAX_CHARGE];
extern float fMassH2O[MAX_CHARGE];
extern float fMassCO[MAX_CHARGE];
    
extern char  cModSites[('Z' - 'A')+1];
extern float fModValues[('Z' - 'A')+1];

extern std::vector<nlValue> nlValuesNterm;
extern std::vector<nlValue> nlValuesCterm;
extern std::vector<nlValue> nlValuesAA;

extern int   *host_iPeakCounts;
extern long  *host_lPeakIndices;
extern std::vector<int>   host_iPeakBins;
extern std::vector<float> host_fPeakInts;
extern cObj  *host_cCandidates;
extern std::stack<cl_mem> unusedBuffers;
extern std::map<long, cl_mem> spectrum2buffer;

extern cl_mem cl_iPeakCounts;
extern cl_mem cl_lPeakIndices;
extern cl_mem cl_iPeakBins;
extern cl_mem cl_fPeakInts;
extern cl_mem cl_fSpectra;
extern cl_mem cl_fScratch;
extern cl_mem cl_bTheoData;
extern cl_mem cl_init_bTheoData;
extern cl_mem cl_init_fSpectra;
extern cl_mem cl_init_fScratch;
extern cl_mem cl_cCandidates;
extern cl_mem cl_fScores;
extern cl_mem cl_mPSMs;
extern cl_mem cl_fNextScores;
extern cl_mem cl_nlValuesNterm;
extern cl_mem cl_nlValuesCterm;
extern cl_mem cl_nlValuesAA;

extern cl_context clContext;
extern cl_platform_id clPlatformID;
extern cl_device_id clDeviceID;
extern cl_program clProgram;
extern cl_command_queue clCommandQueue;
extern long totalScoreTime;
extern long totalReduceTime;
extern long totalBuildTime;
extern long totalMemsetTime;
extern long totalSendTime;
extern long totalTransformTime;
extern long buildLaunches;
extern long scoreKernelLaunches;

//==================================================================================================
//  Functions
//==================================================================================================

// main.c
extern void tempest_exit(int );

// globals.c
extern void initialize_globals(void);
extern void create_kernels(void);
extern void setup_globals(void);
extern void cleanup_globals(void);

// input.c
extern void parse_input( int, char **);

// output.c
extern void write_psms(void);
extern void write_log(void);

// params.c
extern void parse_params();

// experimental.c
extern void collect_msms_spectra(void);
extern void collect_dta_spectra(void);

// theoretical.c
extern void set_residue_masses(void);
extern void search_fasta_database(void);
extern void parse_protein(long, int, char *);

// device.c
extern void initialize_device(void);
extern void setup_device(void);
extern cl_program build_program(cl_context, cl_device_id, const char*);
extern void cleanup_device(void);

// kernels.cu
extern int cl_memset(cl_mem buffer, int c, unsigned long n, cl_event* evt); 
extern void gpu_score_candidates(eObj *);
extern void prebuild_parallel();
extern void prebuild_sequential();
extern void setup_constant_memory(void);

// util.c
extern int n_choose_k(int,int);
extern long mround(long, int);
extern const char *byte_to_binary(unsigned int);
extern int count_set_bits(unsigned int);
extern char* modpeptide(const char*, unsigned int);
extern char* modnpeptide(const char*, int, unsigned int);
extern const char* get_setup_mode(enum SETUP_MODE);
extern const char* get_score_mode(enum SCORE_MODE);
extern unsigned long hash(char*);
extern unsigned long hashn(char*, int);

extern void check_cl_error(char*, int, int, char*);
extern void print_device_info();
extern const char* get_error_string(int);
extern char* strdup_s(char*);

#ifdef _WIN32
extern float roundf(float);
#endif

#endif /*FMACO_H*/
