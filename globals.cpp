#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


//==================================================================================================
//	Variables
//==================================================================================================

struct ARGS   Tempest::args = {0};
struct PARAMS Tempest::params = {0};
struct CONFIG Tempest::config = {0};
struct DATA   Tempest::data = {0};

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
    Tempest::params.fragmentBinOffset = 0.0;
    Tempest::params.bRemovePrecursors = 0;
    Tempest::params.fRemoveMzRanges = 0;
    Tempest::params.iNumRemoveMzRanges = 0;
    Tempest::params.numNormRegions = 10;
    Tempest::params.intensityThreshold = 0.0;
    Tempest::params.xcorrTransformWidth = 0;
    Tempest::params.flankingIntensity = 0.0;
    Tempest::params.numInternalPSMs = 10;
    Tempest::params.numOutputPSMs = 1;
    Tempest::params.bFixDeltaScore = 0;

    //config
    Tempest::config.iPlatform = -1;
    Tempest::config.iDevices = std::vector<unsigned int>();
    Tempest::config.minWorkSize = 32;
    Tempest::config.minScan = 0;
    Tempest::config.maxScan = INT_MAX;
    Tempest::config.parallelReduce = 1;
    Tempest::config.profile = 0;

    //info/status
    Tempest::data.iNumSpectra = 0;
    Tempest::data.lNumMS2Peaks = 0;
    Tempest::data.iNumProteins = 0;
    Tempest::data.lNumPeptides = 0;
    Tempest::data.lNumPSMs = 0;
    Tempest::data.tStart = time(0);
    
    for (int offset=0; ('Z'+offset)<256; offset+=32) {
        // amino-acids monoisotopic mass
        Tempest::params.dMassAA['A'+offset] =  71.03711378;
        Tempest::params.dMassAA['B'+offset] = 114.53493523;
        Tempest::params.dMassAA['C'+offset] = 103.00918451;
        Tempest::params.dMassAA['D'+offset] = 115.02694302;
        Tempest::params.dMassAA['E'+offset] = 129.04259308;
        Tempest::params.dMassAA['F'+offset] = 147.06841390;
        Tempest::params.dMassAA['G'+offset] =  57.02146372;
        Tempest::params.dMassAA['H'+offset] = 137.05891186;
        Tempest::params.dMassAA['I'+offset] = 113.08406396;
        Tempest::params.dMassAA['J'+offset] =   0.0;
        Tempest::params.dMassAA['K'+offset] = 128.09496300;
        Tempest::params.dMassAA['L'+offset] = 113.08406396;
        Tempest::params.dMassAA['M'+offset] = 131.04048463;
        Tempest::params.dMassAA['N'+offset] = 114.04292744;
        Tempest::params.dMassAA['O'+offset] = 114.07931294;
        Tempest::params.dMassAA['P'+offset] =  97.05276384;
        Tempest::params.dMassAA['Q'+offset] = 128.05857750;
        Tempest::params.dMassAA['R'+offset] = 156.10111102;
        Tempest::params.dMassAA['S'+offset] =  87.03202840;
        Tempest::params.dMassAA['T'+offset] = 101.04767846;
        Tempest::params.dMassAA['U'+offset] = 150.04344;
        Tempest::params.dMassAA['V'+offset] =  99.06841390;
        Tempest::params.dMassAA['W'+offset] = 186.07931294;
        Tempest::params.dMassAA['Y'+offset] = 163.06332852;
        Tempest::params.dMassAA['Z'+offset] = 128.55058529;

        Tempest::params.unModAA['A'+offset] = 'A';
        Tempest::params.unModAA['B'+offset] = 'B';
        Tempest::params.unModAA['C'+offset] = 'C';
        Tempest::params.unModAA['D'+offset] = 'D';
        Tempest::params.unModAA['E'+offset] = 'E';
        Tempest::params.unModAA['F'+offset] = 'F';
        Tempest::params.unModAA['G'+offset] = 'G';
        Tempest::params.unModAA['H'+offset] = 'H';
        Tempest::params.unModAA['I'+offset] = 'I';
        Tempest::params.unModAA['J'+offset] = 'J';
        Tempest::params.unModAA['K'+offset] = 'K';
        Tempest::params.unModAA['L'+offset] = 'L';
        Tempest::params.unModAA['M'+offset] = 'M';
        Tempest::params.unModAA['N'+offset] = 'N';
        Tempest::params.unModAA['O'+offset] = 'O';
        Tempest::params.unModAA['P'+offset] = 'P';
        Tempest::params.unModAA['Q'+offset] = 'Q';
        Tempest::params.unModAA['R'+offset] = 'R';
        Tempest::params.unModAA['S'+offset] = 'S';
        Tempest::params.unModAA['T'+offset] = 'T';
        Tempest::params.unModAA['U'+offset] = 'U';
        Tempest::params.unModAA['V'+offset] = 'V';
        Tempest::params.unModAA['W'+offset] = 'W';
        Tempest::params.unModAA['Y'+offset] = 'Y';
        Tempest::params.unModAA['Z'+offset] = 'Z';
    }
    Tempest::params.dMassAA['\0'] = 0.0;
    Tempest::params.unModAA['\0'] = '\0';
}


/*
 *
 */

extern void Tempest::setup_globals()
{
	int i;
	int iNumVariableMods;
	bool bSTPhos;

	// MS/MS scans
	if (Tempest::params.bPrecursorTolerancePPM) Tempest::data.fPrecursorBinWidth = Tempest::params.maxPeptideMass * Tempest::params.fPrecursorTolerance / 1000000.0f;
	else Tempest::data.fPrecursorBinWidth = Tempest::params.fPrecursorTolerance;

	Tempest::data.iNumPrecursorBins = (int) ceil(Tempest::params.maxPeptideMass / Tempest::data.fPrecursorBinWidth);

	// MS/MS Scan Index
	if (0 == (Tempest::data.eScanIndex = (eObj**) calloc(Tempest::data.iNumPrecursorBins, sizeof(eObj)))) {
		fprintf(stderr, "\nERROR! Unable to malloc %d bytes for MS/MS Scan Index\n", (int) (Tempest::data.iNumPrecursorBins * sizeof(eObj)));
		Tempest::tempest_exit(EXIT_FAILURE);
	}

	// MS/MS spectrum
    // currently cross correlation window is # of bins and does not scale to bin width
	//Tempest::data.iCrossCorrelationWidth = (int) CROSS_CORRELATION_WINDOW / Tempest::params.fFragmentTolerance;
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
	// 	for (i=0; i<Tempest::data.iNumPrecursorBins;i++) {
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
