#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


//==================================================================================================
//	Variables
//==================================================================================================
namespace Tempest {
    struct ARGS args;
    struct PARAMS params;
    struct CONFIG config;
    struct INFO tempest;
    struct GPUINFO gpu_info;

    // MS/MS scans
    eObj **eScanIndex;
    bool **bScanIndex;

    // proteins
    char *sProteinReferences;

    // masses
    double dMassAA[256] = {0.0};

    // mods
    int   numAAModSites[256] = {0};
    char  cModSites[256][5] = {'\0'};
    float fModValues[256][5] = {0.0};
    char  unModAA[256];
    int   numNtermModSites;
    char  ntermModSymbols[5] = {'\0'};
    float ntermModMasses[5] = {0};
    int   numCtermModSites = 0;
    char  ctermModSymbols[5] = {'\0'};
    float ctermModMasses[5] = {0};

    // neutral losses
    //for sorting neutral_loss parameter lines
    std::vector<nlValue> nlValuesNterm;
    std::vector<nlValue> nlValuesCterm;
    std::vector<nlValue> nlValuesAA;

    // host memory
    int   *host_iPeakCounts;
    long  *host_lPeakIndices;
    cObj  *host_cCandidates;

    bool* bDigestSites;
    bool* bDigestNoSites;
}

//==================================================================================================
//	Functions
//==================================================================================================

/*
 * Initialize global variables
 */

extern void Tempest::initialize_globals()
{
	int i;

    //command line arguments
    Tempest::args.sSpectra = 0;
    Tempest::args.sFasta = 0;
    Tempest::args.sParams = 0;
    Tempest::args.sConfig = 0;
    Tempest::args.sOut = 0;
    Tempest::args.bSplitOutput = 0;
    Tempest::args.iPrintLevel = 1;
    Tempest::args.bPrintCandidates = 0;
    Tempest::args.bNoDigest = 0;
    Tempest::args.bNoMatch = 0;
    Tempest::args.bNoScore = 0;
    Tempest::args.deviceOverride = 0;
    Tempest::args.hostMemOverride = 0;
    Tempest::args.deviceMemOverride = 0;
	
    //params
    Tempest::params.sDigestSites = 0;
    Tempest::params.sDigestNoSites = 0;
    Tempest::params.iDigestSpecificity = -1;
    Tempest::params.iDigestMaxMissedCleavages = 2;
    Tempest::params.iMinPeptideLength = 5;
    Tempest::params.iMaxPeptideLength = 64;
    Tempest::params.minPeptideMass = 400.0;
    Tempest::params.maxPeptideMass = 5000.0;
    Tempest::params.tMods = 0;
    Tempest::params.iNumMods = 0;
    Tempest::params.iModificationsMax = 3;
    Tempest::params.numAANL = 0;
    Tempest::params.numNtermNL = 0;
    Tempest::params.numCtermNL = 0;
    Tempest::params.useAIons = false;
    Tempest::params.useBIons = false;
    Tempest::params.useCIons = false;
    Tempest::params.useXIons = false;
    Tempest::params.useYIons = false;
    Tempest::params.useZIons = false;
    Tempest::params.dPeptideNtermMass = 0.0;
    Tempest::params.dPeptideCtermMass = 0.0;
    Tempest::params.dProteinNtermMass = 0.0;
    Tempest::params.dProteinCtermMass = 0.0;
    Tempest::params.fPrecursorTolerance = -1;
    Tempest::params.fFragmentTolerance = -1;
    Tempest::params.fragmentBinOffset = 0.5;
    Tempest::params.bRemovePrecursors = 0;
    Tempest::params.fRemoveMzRanges = 0;
    Tempest::params.iNumRemoveMzRanges = 0;
    Tempest::params.numNormRegions = 0;
    Tempest::params.intensityThreshold = 0.0;
    Tempest::params.xcorrTransformWidth = 0;
    Tempest::params.flankingIntensity = 0.0;
    Tempest::params.numInternalPSMs = 10;
    Tempest::params.numOutputPSMs = 1;
    Tempest::params.bFixDeltaScore = 0;

    //config
    Tempest::config.iPlatform = -1;
    Tempest::config.iDevices = std::vector<unsigned int>();
    Tempest::config.minScan = 0;
    Tempest::config.maxScan = INT_MAX;

    //info/status
    Tempest::tempest.iNumSpectra = 0;
    Tempest::tempest.lNumMS2Peaks = 0;
    Tempest::tempest.iNumProteins = 0;
    Tempest::tempest.lNumPeptides = 0;
    Tempest::tempest.lNumPSMs = 0;
    Tempest::tempest.tStart = time(0);

    // GPU profiling
    Tempest::gpu_info.bTimeGPU = 0;
    Tempest::gpu_info.bTimeSend = 0;
    Tempest::gpu_info.bTimeKernel = 0;
    Tempest::gpu_info.bTimeIdle = 0;
    Tempest::gpu_info.bTimePrebuild = 0;
    Tempest::gpu_info.bCountSyncBlocks = 0;
	
    Tempest::gpu_info.fTotalTimeGPU = 0.0f;
    Tempest::gpu_info.fTotalTimeSend = 0.0f;
    Tempest::gpu_info.fTotalTimeKernel = 0.0f;
    Tempest::gpu_info.fTotalTimePrebuild = 0.0f;
    Tempest::gpu_info.fTotalTimeIdle = 0.0f;
    Tempest::gpu_info.iNumSyncBlocks = 0;
    Tempest::gpu_info.iNumScoringKernels = 0;
    
    for (int offset=0; ('Z'+offset)<256; offset+=32) {
        // amino-acids monoisotopic mass
        Tempest::dMassAA['A'+offset] =  71.03711378;
        Tempest::dMassAA['B'+offset] = 114.53493523;
        Tempest::dMassAA['C'+offset] = 103.00918451;
        Tempest::dMassAA['D'+offset] = 115.02694302;
        Tempest::dMassAA['E'+offset] = 129.04259308;
        Tempest::dMassAA['F'+offset] = 147.06841390;
        Tempest::dMassAA['G'+offset] =  57.02146372;
        Tempest::dMassAA['H'+offset] = 137.05891186;
        Tempest::dMassAA['I'+offset] = 113.08406396;
        Tempest::dMassAA['J'+offset] =   0.0;
        Tempest::dMassAA['K'+offset] = 128.09496300;
        Tempest::dMassAA['L'+offset] = 113.08406396;
        Tempest::dMassAA['M'+offset] = 131.04048463;
        Tempest::dMassAA['N'+offset] = 114.04292744;
        Tempest::dMassAA['O'+offset] = 114.07931294;
        Tempest::dMassAA['P'+offset] =  97.05276384;
        Tempest::dMassAA['Q'+offset] = 128.05857750;
        Tempest::dMassAA['R'+offset] = 156.10111102;
        Tempest::dMassAA['S'+offset] =  87.03202840;
        Tempest::dMassAA['T'+offset] = 101.04767846;
        Tempest::dMassAA['U'+offset] = 150.04344;
        Tempest::dMassAA['V'+offset] =  99.06841390;
        Tempest::dMassAA['W'+offset] = 186.07931294;
        Tempest::dMassAA['Y'+offset] = 163.06332852;
        Tempest::dMassAA['Z'+offset] = 128.55058529;

        Tempest::unModAA['A'+offset] = 'A';
        Tempest::unModAA['B'+offset] = 'B';
        Tempest::unModAA['C'+offset] = 'C';
        Tempest::unModAA['D'+offset] = 'D';
        Tempest::unModAA['E'+offset] = 'E';
        Tempest::unModAA['F'+offset] = 'F';
        Tempest::unModAA['G'+offset] = 'G';
        Tempest::unModAA['H'+offset] = 'H';
        Tempest::unModAA['I'+offset] = 'I';
        Tempest::unModAA['J'+offset] = 'J';
        Tempest::unModAA['K'+offset] = 'K';
        Tempest::unModAA['L'+offset] = 'L';
        Tempest::unModAA['M'+offset] = 'M';
        Tempest::unModAA['N'+offset] = 'N';
        Tempest::unModAA['O'+offset] = 'O';
        Tempest::unModAA['P'+offset] = 'P';
        Tempest::unModAA['Q'+offset] = 'Q';
        Tempest::unModAA['R'+offset] = 'R';
        Tempest::unModAA['S'+offset] = 'S';
        Tempest::unModAA['T'+offset] = 'T';
        Tempest::unModAA['U'+offset] = 'U';
        Tempest::unModAA['V'+offset] = 'V';
        Tempest::unModAA['W'+offset] = 'W';
        Tempest::unModAA['Y'+offset] = 'Y';
        Tempest::unModAA['Z'+offset] = 'Z';
    }
    Tempest::dMassAA['\0'] = 0.0;
    Tempest::unModAA['\0'] = '\0';
}


/*
 *
 */

extern void Tempest::setup_globals()
{
	int i;
	int iNumVariableMods;
	bool bSTPhos;

    //MEA: took this out
	// count variable mods and look for Serine/Threonine phosphorylation mod
	// iNumVariableMods = 0;
	// bSTPhos = 0;
	// for (i=0; i<Tempest::params.iNumMods; i++) {
	// 	if (Tempest::params.tMods[i].cSymbol) {
	// 		iNumVariableMods++;
	// 		if (abs(Tempest::params.tMods[i].dMassDiff-80) < 1.0 && (Tempest::params.tMods[i].cAminoAcid == 'S' || Tempest::params.tMods[i].cAminoAcid == 'T'))
	// 			bSTPhos = 1;
	// 	}
	// }

	//constrain max modifications
	// if (iNumVariableMods == 0) {
	// 	Tempest::params.iModificationsMax = 0;
	// 	if (Tempest::params.cVariableNtermSymbol) Tempest::params.iModificationsMax += 1;
	// 	if (Tempest::params.cVariableCtermSymbol) Tempest::params.iModificationsMax += 1;
	// }

	// MS/MS scans
	if (Tempest::params.bPrecursorTolerancePPM) Tempest::tempest.fPrecursorBinWidth = Tempest::params.maxPeptideMass * Tempest::params.fPrecursorTolerance / 1000000.0f;
	else Tempest::tempest.fPrecursorBinWidth = Tempest::params.fPrecursorTolerance;

	Tempest::tempest.iNumPrecursorBins = (int) ceil(Tempest::params.maxPeptideMass / Tempest::tempest.fPrecursorBinWidth);

	// MS/MS Scan Index
	if (0 == (Tempest::eScanIndex = (eObj**) calloc(Tempest::tempest.iNumPrecursorBins, sizeof(eObj)))) {
		fprintf(stderr, "\nERROR! Unable to malloc %d bytes for MS/MS Scan Index\n", (int) (Tempest::tempest.iNumPrecursorBins * sizeof(eObj)));
		Tempest::tempest_exit(EXIT_FAILURE);
	}
	
	//MS/MS Scan Index (booleans for mods)
	if (0 == (Tempest::bScanIndex = (bool**) malloc((Tempest::params.iModificationsMax+1) * sizeof(bool*)))) {
		fprintf(stderr, "\nERROR! Unable to allocate %d bytes for MS/MS Scan Modification Index\n", (int) ((Tempest::params.iModificationsMax+1) * sizeof(bool*)));
		Tempest::tempest_exit(EXIT_FAILURE);
	}
    
	for (i=0; i<=Tempest::params.iModificationsMax; i++) {
		if (0 == (Tempest::bScanIndex[i] = (bool*) calloc(Tempest::tempest.iNumPrecursorBins, sizeof(bool)))) {
			fprintf(stderr, "\nERROR! Unable to allocate %d bytes for MS/MS Scan Index w/ %d mods\n", (int) (Tempest::tempest.iNumPrecursorBins * sizeof(bool)), i);
			Tempest::tempest_exit(EXIT_FAILURE);
		}
	}

	// MS/MS spectrum
	//Tempest::tempest.iCrossCorrelationWidth = (int) CROSS_CORRELATION_WINDOW / Tempest::params.fFragmentTolerance;

    //setup enzyme sites
    Tempest::bDigestSites = (bool*) calloc(256, sizeof(bool));
    Tempest::bDigestNoSites = (bool*) calloc(256, sizeof(bool));
        
    for (unsigned int i=0; i<strlen(Tempest::params.sDigestSites); i++) {
        Tempest::bDigestSites[Tempest::params.sDigestSites[i]] = 1;
    }

    for (unsigned int i=0; i<strlen(Tempest::params.sDigestNoSites); i++) {
        Tempest::bDigestNoSites[Tempest::params.sDigestNoSites[i]] = 1;
    }
}

/*
 * Free memory
 */
/*
extern void cleanup_globals() {
	int i;
	eObj *e, *eNext;

	// args
	//safe_free(Tempest::args.sSpectra);
	//safe_free(Tempest::args.sFasta);
	//safe_free(Tempest::args.sParams);
	//safe_free(Tempest::args.sOut);

	// params
	//safe_free(Tempest::params.tMods);
	Tempest::params.iNumMods=0;

	//safe_free(Tempest::params.fRemoveMzRanges);
	Tempest::params.iNumRemoveMzRanges=0;

	// // scans
	// if (eScanIndex) {
	// 	for (i=0; i<Tempest::tempest.iNumPrecursorBins;i++) {
	// 		e = eScanIndex[i];
	// 		while (e) {
	// 			safe_free(e->sName);
	// 			//cudaEventDestroy(e->cudaEventSent);
    //             clReleaseEvent(e->clEventSent);
	// 			eNext = e->next;
	// 			free(e);
	// 			e = eNext;
	// 		}
	// 	}

	// 	free(eScanIndex);
	// 	eScanIndex = 0;
	// }

	// // scan Index
	// if (bScanIndex) {
	// 	for (i=0; i<=Tempest::params.iModificationsMax; i++)
	// 		safe_free(bScanIndex[i]);

	// 	free(bScanIndex);
	// 	bScanIndex = 0;
	// }

	// proteins
    Tempest::safe_free(sProteinReferences);

    //MEA
	// candidates (pinned host memory)
	//safe_cudaFreeHost(host_cCandidates);
}
*/
