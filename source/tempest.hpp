#ifndef TEMPEST_HPP
#define TEMPEST_HPP

#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
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
#include <set>
#include <map>

typedef off_t f_off; //for MSToolkit
#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"

#include "constants.h" //molecular masses

#define VERSION_STRING "2.0"
#define LICENSE_STRING ""
#define AUTHORS_STRING "Mark Adamo, Jeffrey Milloy, Brendan Faherty"

//==================================================================================================
//  Settings and Defaults
//==================================================================================================

#define STRING_SIZE 4096
#define CASE_SHIFT 32

#define KB 1024
#define MB 1048576
#define GB 1073741824

// params
#define DEFAULT_PARAMS_FILE "tempest.method"
#define DEFAULT_CONFIG_FILE "tempest.config"
// kernels
#define KERNEL_FILE "kernels.cl"

// digestion
#define MAX_LENGTH_REFERENCE 64
#define MIN_PEPTIDE_LENGTH 1
#define MAX_PEPTIDE_LENGTH 64
#define MAX_MOD_TARGETS_PER_PEPTIDE 64
#define INITIAL_PROTEIN_BUFFER 2048

#define INT_MAX 2147483647;

#define MIN_SCORE -1000.0f

//==================================================================================================
// Macros
//==================================================================================================

#define safe_free(p) if (p) { free(p); p=0; }

//==================================================================================================
//  Structures
//==================================================================================================

class Device;

// Candidate Peptide
typedef struct cobject {
    int    iProtein;
    float  fPeptideMass;

    char   decoy; // bool (for internal decoy_search)
    unsigned char iPeptideLength;
    char   ntermMod; //values from 0 to 5
    char   ctermMod; //values from 0 to 5
    char   cBefore;
    char   cAfter;
    
    unsigned char sPeptide[MAX_PEPTIDE_LENGTH];
} cObj;

// Peptide-Spectrum Match
typedef struct mobject {
    int    iPeptideLength;
    int    iNumOccurrences;
    float  fPeptideMass;
    float  fScore;
    int    iProtein;

    char   decoy; // bool (for internal decoy_search)
    char   ntermMod; //values from 0 to 5
    char   ctermMod; //values from 0 to 5
    char   cBefore;
    char   cAfter;

    unsigned char sPeptide[MAX_PEPTIDE_LENGTH];
} mObj;

// Observed Scan
typedef struct eobject {
    long     lIndex;
    char*    sName;
    double   dPrecursorMass;
    int      iPrecursorCharge;
    cObj*    candidateBuffer;
    size_t   candidateBufferSize;
    unsigned int iNumBufferedCandidates;
    unsigned int iNumCandidates;
    cl_event clEventSent;
    Device* device;
    //cudaEvent_t cudaEventSent;
    
    struct eobject* next;
} eObj;

//Mass delta and score weighting for a neutral loss event
typedef struct nlValue {
    int   numSites;
    char  hasAA[256]; //boolean
    char  modNum; //values from 0 to 5
    float massDelta;
    float weighting;
} nlValue;

typedef struct MOD {
    unsigned char    cAminoAcid;
    int     modInd;
    char    cSymbol;
    double  dMassDiff;
} mod_t;

struct ARGS {
    char *sFasta;
    char *sSpectra;
    char *sOut;
    const char *sParams;
    const char *sConfig;

    int iPrintLevel;
    bool bPrintCandidates;

    bool deviceOverride;
    bool hostMemOverride;
    bool deviceMemOverride;
};

struct PARAMS {
    //files
    char *sFasta;
    char *sSpectra;
    char *sOut;

    //decoys
    bool decoySearch;
    
    //static AA masses
    double dMassAA[256];
    
    // digestion
    char sDigestSites[26];
    char sDigestNoSites[26];
    bool bDigestSites[256];
    bool bDigestNoSites[256];
    int iDigestOffset;
    int iDigestSpecificity;
    int iDigestMaxMissedCleavages;
    int iMinPeptideLength;
    int iMaxPeptideLength;
    float minPeptideMass;
    float maxPeptideMass;
    
    // modifications
    mod_t *tMods;
    int    iNumMods;
    int    iModificationsMax;

    //neutral loss counts
    int   numAANL;
    int   numNtermNL;
    int   numCtermNL;
    
    //neutral losses
    int numAAModSites[256];
    char cModSites[256][5];
    float fModValues[256][5];
    char unModAA[256];
    int numNtermModSites;
    char ntermModSymbols[5];
    float ntermModMasses[5];
    bool ntermModProtein[5];
    int numCtermModSites;
    char ctermModSymbols[5];
    float ctermModMasses[5];
    bool ctermModProtein[5];
    std::vector<nlValue> nlValuesNterm;
    std::vector<nlValue> nlValuesCterm;
    std::vector<nlValue> nlValuesAA;
    
    //ion series
    float AIons;
    float BIons;
    float CIons;
    float XIons;
    float YIons;
    float ZIons;

    //MS level
    std::set<int> msLevels;

    double dPeptideNtermMass;
    double dPeptideCtermMass;
    double dProteinNtermMass;
    double dProteinCtermMass;

    // precursor (MS1) tolerance
    float  fPrecursorTolerance;
    bool   bPrecursorTolerancePPM;

    // fragment (MS2) tolerance
    float  fFragmentTolerance;
    bool   bFragmentTolerancePPM;
    float  fragmentBinOffset;

    // spectra processing
    bool   bRemovePrecursors;
    float *fRemoveMzRanges;
    int    iNumRemoveMzRanges;
    int    numNormRegions;
    float  intensityThreshold;

    // similarity scoring
    int   xcorrTransformWidth;
    float flankingIntensity;
    int   numInternalPSMs;

    // output
    int   numOutputPSMs;
    bool  bFixDeltaScore;
};

struct CONFIG {
    unsigned int iPlatform;
    std::vector<unsigned int> iDevices;
    // scan range
    int minScan;
    int maxScan;
    size_t maxHostMem;
    size_t maxDeviceMem;
    int minWorkSize;
    bool parallelReduce;
    bool profile;
};

struct DATA {
    // MS/MS scan index
    eObj **eScanIndex;
    std::vector<eObj*> eScans;
    
    // MS/MS scan info
    int    iNumSpectra;
    float  fPrecursorBinWidth;
    int    iNumPrecursorBins;
    float  fMinPrecursorMass;
    float  fMaxPrecursorMass;

    // MS/MS peak data
    int   *host_iPeakCounts;
    long  *host_lPeakIndices;
    cObj  *host_cCandidates;
    std::vector<int>   host_iPeakBins;
    std::vector<float> host_fPeakInts;

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
    char   *sProteinReferences;

    // Runtime info
    time_t tStart;
};

class Device {
private:
    int platformID;
    int deviceID;
    int deviceInd; //index in devices vector
    unsigned int minScan;
    unsigned int maxScan;
    
    cl_mem cl_iPeakCounts;
    cl_mem cl_lPeakIndices;
    cl_mem cl_iPeakBins;
    cl_mem cl_fPeakInts;
    cl_mem cl_fSpectra;
    cl_mem cl_init_fSpectra;
    cl_mem cl_cCandidates;
    cl_mem cl_fScores;
    cl_mem cl_mPSMs;
    cl_mem cl_nlValuesNterm;
    cl_mem cl_nlValuesCterm;
    cl_mem cl_nlValuesAA;
    std::stack<cl_mem> unusedBuffers;

    // global kernels
    cl_kernel __cl_build;
    cl_kernel __cl_transform;
    cl_kernel __cl_score;
    cl_kernel __cl_reduce_scores;
    cl_kernel __cl_memset;
    //kernel local work dimensions
    size_t build_size;
    size_t transform_size;
    size_t score_size;
    size_t reduce_scores_multiple;
    size_t reduce_scores_size_local;
    size_t reduce_scores_size;
    size_t reduce_scores_size_max;
    size_t memset_size;
    //buffer size
    size_t candidateBufferSize;

    std::map<long, cl_mem> spectrum2buffer;

    cl_event scoreEvent;
    cl_event reduceEvent;
    cl_event buildEvent;
    cl_event memsetEvent;
    cl_event transformEvent;

    //profiling
    long totalScoreTime;
    long totalReduceTime;
    long totalBuildTime;
    long totalTransformTime;
    long totalMemsetTime;
    long buildLaunches;
    long scoreKernelLaunches;
    long lastBuildIndex;

    cl_context clContext;
    cl_platform_id clPlatformID;
    cl_device_id clDeviceID;
    cl_program clProgram;
    cl_kernel clKernel;
    cl_command_queue clCommandQueue;

    cl_uint  iMaxComputeUnits;
    cl_uint  iAddressBits;
    cl_bool  bUnifiedMemory;
    cl_ulong lMaxWorkGroupSize;
    cl_ulong lGlobalMemSize;
    cl_ulong lLocalMemSize;
    cl_ulong lMaxMemAllocSize;

    int cl_memset(cl_mem buffer, int c, unsigned long n, cl_event* evt);
    void setup_constant_memory();
    void create_kernels();

public:
    Device(int platformID, int deviceInd);
    void setup(unsigned int minScan, unsigned int maxScan);
    void scoreCandidates(eObj *e);
    void finish();
    cl_event newEvent();
    int get_mPSMs(mObj* destination);
    int get_fNextScores(float* destination);
    long totalSendTime;
    void printProfilingData();
};

//==================================================================================================
//  Global namespace variables
//==================================================================================================

namespace Tempest {
    extern std::vector<Device*> devices;

    extern struct ARGS args;
    extern struct PARAMS params;
    extern struct CONFIG config;
    extern struct DATA data;

    //==================================================================================================
    //  Functions
    //==================================================================================================

    // main.cpp
    extern void tempest_exit(int);

    // globals.cpp
    extern void initialize_globals(void);
    extern void create_kernels(void);
    extern void setup_globals(void);
    extern void cleanup_globals(void);

    // input.cpp
    extern void parse_input(int, char **);
    extern void usage();

    // output.cpp
    extern void write_psms(void);
    extern void write_log(int, char**);
    // params.cpp
    extern void parse_params();

    // config.cpp
    extern void parse_config();

    // experimental.cpp
    extern void collect_msms_spectra(void);
    extern void collect_dta_spectra(void);

    // theoretical.cpp
    extern void set_residue_masses(void);
    extern void search_fasta_database(void);
    extern void parse_protein(long, int, char *);

    // device.cpp
    extern void initialize_device(void);
    extern void setup_device(void);
    extern cl_program build_program(cl_context, cl_device_id, const char*);
    extern void cleanup_device(void);

    // util.cpp
    extern long mround(long, int);
    extern unsigned char toMod(char c, int modInd);
    extern int getModInd(unsigned char c);
    extern bool backboneMatch(mObj m1, mObj m2);

    extern void check_cl_error(const char*, int, int, const char*);
    extern void print_device_info();
    extern const char* get_error_string(int);
    extern char* strdup_s(char*);
}

#endif
