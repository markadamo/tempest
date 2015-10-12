#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


//==================================================================================================
//	Variables
//==================================================================================================

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
double dMassAA[('z' - 'A')+1];
float fMassProton[MAX_CHARGE];
float fMassNH3[MAX_CHARGE];
float fMassH2O[MAX_CHARGE];
float fMassCO[MAX_CHARGE];

// mods
char  cModSites[('Z' - 'A')+1];
float fModValues[('Z' - 'A')+1];

// neutral losses 
//for sorting neutral_loss parameter lines
std::vector<nlValue> nlValuesNterm;
std::vector<nlValue> nlValuesCterm;
std::vector<nlValue> nlValuesAA;

// host memory
int   *host_iPeakCounts;
long  *host_lPeakIndices;
cObj  *host_cCandidates;

//config constants
size_t MAX_BLOCKDIM = 32;
size_t DEFAULT_CANDIDATE_BUFFER_SIZE = 1024;
size_t DEFAULT_BLOCKDIM_SCORE = 32;
size_t BLOCKDIM_PREBUILD = 1;
size_t BLOCKDIM_PREBUILD_TRANSFORMATION = 1;
size_t BLOCKDIM_BUILD = 32;
size_t BLOCKDIM_TRANSFORM = 32;
size_t BLOCKDIM_REDUCE = 1024;
size_t CL_ONE = 256;

//==================================================================================================
//	Functions
//==================================================================================================

/*
 * Initialize global variables
 */

extern void initialize_globals()
{
	int i;
	
	//command line arguments
	args.sSpectra = 0;
	args.sFasta = 0;
	args.sParams = 0;
	args.sOut = 0;
	args.bSplitOutput = 0;
	args.iPrintLevel = 1;
	args.bPrintCandidates = 0;
	args.bNoDigest = 0;
	args.bNoMatch = 0;
	args.bNoScore = 0;
	
	//params
	params.sDigestSites = 0;
	params.sDigestNoSites = 0;
	params.iDigestSpecificity = -1;
	params.iDigestMaxMissedCleavages = 2;
	params.iMinPeptideLength = MIN_PEPTIDE_LENGTH;
	params.iMaxPeptideLength = MAX_PEPTIDE_LENGTH;
	params.tMods = 0;
	params.iNumMods = 0;
	params.iModificationsMax = 3;
    params.numAANL = 0;
    params.numNtermNL = 0;
    params.numCtermNL = 0;
	params.dPeptideNtermMass = 0.0;
	params.dPeptideCtermMass = 0.0;
	params.dProteinNtermMass = 0.0;
	params.dProteinCtermMass = 0.0;
	params.dVariableNtermMassDiff = 0.0;
	params.dVariableCtermMassDiff = 0.0;
	params.cVariableNtermSymbol = '\0';
	params.cVariableCtermSymbol = '\0';
	params.fPrecursorTolerance = -1;
	params.fFragmentTolerance = -1;
	params.bRemovePrecursors = 0;
	params.fRemoveMzRanges = 0;
	params.iNumRemoveMzRanges = 0;
	params.bCrossCorrelation = 0;
	params.bTheoreticalFlanking = 1;
	params.iNumOutputPSMs = 1;
	params.bFixDeltaScore = 0;

	//config
    config.iPlatform = 0;
	config.iDevices = std::vector<unsigned int>();
	config.iCandidateBufferSize = DEFAULT_CANDIDATE_BUFFER_SIZE;
	config.iScoreBlockDim = DEFAULT_BLOCKDIM_SCORE;
    config.minScan = 0;
    config.maxScan = INT_MAX;

	//info/status
	tempest.iNumSpectra = 0;
	tempest.lNumMS2Peaks = 0;
	tempest.iNumProteins = 0;
	tempest.lNumPeptides = 0;
	tempest.lNumPSMs = 0;
	tempest.tStart = time(0);

	// GPU profiling
	gpu_info.bTimeGPU = 0;
	gpu_info.bTimeSend = 0;
	gpu_info.bTimeKernel = 0;
	gpu_info.bTimeIdle = 0;
	gpu_info.bTimePrebuild = 0;
	gpu_info.bCountSyncBlocks = 0;
	
	gpu_info.fTotalTimeGPU = 0.0f;
	gpu_info.fTotalTimeSend = 0.0f;
	gpu_info.fTotalTimeKernel = 0.0f;
	gpu_info.fTotalTimePrebuild = 0.0f;
	gpu_info.fTotalTimeIdle = 0.0f;
	gpu_info.iNumSyncBlocks = 0;
	gpu_info.iNumScoringKernels = 0;

	// Initialize the AA modifications
	for(i=aa2i('A');i<=aa2i('Z');i++) {
		cModSites[i] = '\0';
		fModValues[i] = 0.0;
	}
	
	// initialize masses
	for(i=0;i<MAX_CHARGE;i++)
	{
		fMassProton[i] = (float) PROTON_MASS * (float) (i);
		fMassNH3[i] = (float) NH3_MASS / (float) (i+1);
		fMassH2O[i] = (float) H2O_MASS / (float) (i+1);
		fMassCO[i] = (float) CO_MASS / (float) (i+1);
	}

	// amino-acids monoisotopic mass
	dMassAA[aa2i('A')] =  71.03711378;
	dMassAA[aa2i('B')] = 114.53493523;
	dMassAA[aa2i('C')] = 103.00918451;
	dMassAA[aa2i('D')] = 115.02694302;
	dMassAA[aa2i('E')] = 129.04259308;
	dMassAA[aa2i('F')] = 147.06841390;
	dMassAA[aa2i('G')] =  57.02146372;
	dMassAA[aa2i('H')] = 137.05891186;
	dMassAA[aa2i('I')] = 113.08406396;
	dMassAA[aa2i('J')] =   0.0;
	dMassAA[aa2i('K')] = 128.09496300;
	dMassAA[aa2i('L')] = 113.08406396;
	dMassAA[aa2i('M')] = 131.04048463;
	dMassAA[aa2i('N')] = 114.04292744;
	dMassAA[aa2i('O')] = 114.07931294;
	dMassAA[aa2i('P')] =  97.05276384;
	dMassAA[aa2i('Q')] = 128.05857750;
	dMassAA[aa2i('R')] = 156.10111102;
	dMassAA[aa2i('S')] =  87.03202840;
	dMassAA[aa2i('T')] = 101.04767846;
	dMassAA[aa2i('U')] = 150.04344;
	dMassAA[aa2i('V')] =  99.06841390;
	dMassAA[aa2i('W')] = 186.07931294;
	dMassAA[aa2i('Y')] = 163.06332852;
	dMassAA[aa2i('Z')] = 128.55058529;

	// modified mass
	dMassAA[aa2i('a')] =  71.03711378;
	dMassAA[aa2i('b')] = 114.53493523;
	dMassAA[aa2i('c')] = 103.00918451;
	dMassAA[aa2i('d')] = 115.02694302;
	dMassAA[aa2i('e')] = 129.04259308;
	dMassAA[aa2i('f')] = 147.06841390;
	dMassAA[aa2i('g')] =  57.02146372;
	dMassAA[aa2i('h')] = 137.05891186;
	dMassAA[aa2i('i')] = 113.08406396;
	dMassAA[aa2i('j')] =   0.0;
	dMassAA[aa2i('k')] = 128.09496300;
	dMassAA[aa2i('l')] = 113.08406396;
	dMassAA[aa2i('m')] = 131.04048463;
	dMassAA[aa2i('n')] = 114.04292744;
	dMassAA[aa2i('o')] = 114.07931294;
	dMassAA[aa2i('p')] =  97.05276384;
	dMassAA[aa2i('q')] = 128.05857750;
	dMassAA[aa2i('r')] = 156.10111102;
	dMassAA[aa2i('s')] =  87.03202840;
	dMassAA[aa2i('t')] = 101.04767846;
	dMassAA[aa2i('u')] = 150.04344;
	dMassAA[aa2i('v')] =  99.06841390;
	dMassAA[aa2i('w')] = 186.07931294;
	dMassAA[aa2i('y')] = 163.06332852;
	dMassAA[aa2i('z')] = 128.55058529;
}


/*
 *
 */

extern void setup_globals()
{
	int i;
	int iNumVariableMods;
	bool bSTPhos;

	// count variable mods and look for Serine/Threonine phosphorylation mod
	iNumVariableMods = 0;
	bSTPhos = 0;
	for (i=0; i<params.iNumMods; i++) {
		if (params.tMods[i].cSymbol) {
			iNumVariableMods++;
			if (abs(params.tMods[i].dMassDiff-80) < 1.0 && (params.tMods[i].cAminoAcid == 'S' || params.tMods[i].cAminoAcid == 'T'))
				bSTPhos = 1;
		}
	}

	//constrain max modifications
	if (iNumVariableMods == 0) {
		params.iModificationsMax = 0;
		if (params.cVariableNtermSymbol) params.iModificationsMax += 1;
		if (params.cVariableCtermSymbol) params.iModificationsMax += 1;
	}

	// MS/MS scans
	if (params.bPrecursorTolerancePPM) tempest.fPrecursorBinWidth = MAX_DIGEST_MASS * params.fPrecursorTolerance / 1000000.0f;
	else tempest.fPrecursorBinWidth = params.fPrecursorTolerance;

	tempest.iNumPrecursorBins = (int) ceil(MAX_DIGEST_MASS / tempest.fPrecursorBinWidth);

	// MS/MS Scan Index
	if (0 == (eScanIndex = (eObj**) calloc(tempest.iNumPrecursorBins, sizeof(eObj)))) {
		fprintf(stderr, "\nERROR! Unable to malloc %d bytes for MS/MS Scan Index\n", (int) (tempest.iNumPrecursorBins * sizeof(eObj)));
		tempest_exit(EXIT_FAILURE);
	}
	
	//MS/MS Scan Index (booleans for mods)
	if (0 == (bScanIndex = (bool**) malloc((params.iModificationsMax+1) * sizeof(bool*)))) {
		fprintf(stderr, "\nERROR! Unable to allocate %d bytes for MS/MS Scan Modification Index\n", (int) ((params.iModificationsMax+1) * sizeof(bool*)));
		tempest_exit(EXIT_FAILURE);
	}
    
	for (i=0; i<=params.iModificationsMax; i++) {
		if (0 == (bScanIndex[i] = (bool*) calloc(tempest.iNumPrecursorBins, sizeof(bool)))) {
			fprintf(stderr, "\nERROR! Unable to allocate %d bytes for MS/MS Scan Index w/ %d mods\n", (int) (tempest.iNumPrecursorBins * sizeof(bool)), i);
			tempest_exit(EXIT_FAILURE);
		}
	}

	// MS/MS spectrum
	tempest.iCrossCorrelationWidth = (int) CROSS_CORRELATION_WINDOW / params.fFragmentTolerance;

	// config
	config.iScoreNumBlocks = config.iCandidateBufferSize / config.iScoreBlockDim;
}

/*
 * Free memory
 */

extern void cleanup_globals()
{
	int i;
	eObj *e, *eNext;

	// args
	//safe_free(args.sSpectra);
	//safe_free(args.sFasta);
	//safe_free(args.sParams);
	//safe_free(args.sOut);

	// params
	//safe_free(params.tMods);
	params.iNumMods=0;

	//safe_free(params.fRemoveMzRanges);
	params.iNumRemoveMzRanges=0;

	// scans
	if (eScanIndex) {
		for (i=0; i<tempest.iNumPrecursorBins;i++) {
			e = eScanIndex[i];
			while (e) {
				safe_free(e->sName);
				//cudaEventDestroy(e->cudaEventSent);
                clReleaseEvent(e->clEventSent);
				eNext = e->next;
				free(e);
				e = eNext;
			}
		}

		free(eScanIndex);
		eScanIndex = 0;
	}

	// scan Index
	if (bScanIndex) {
		for (i=0; i<=params.iModificationsMax; i++)
			safe_free(bScanIndex[i]);

		free(bScanIndex);
		bScanIndex = 0;
	}

	// proteins
	safe_free(sProteinReferences);

    //MEA
	// candidates (pinned host memory)
	//safe_cudaFreeHost(host_cCandidates);
}
