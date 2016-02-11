#include "tempest.hpp"

#ifdef MAC
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


//==================================================================================================
//	Variables
//==================================================================================================

// zero-initialize variables and arrays in structs
struct ARGS   Tempest::args   = {0};
struct PARAMS Tempest::params = {0};
struct CONFIG Tempest::config = {0};
struct DATA   Tempest::data   = {0};

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
    Tempest::args.iPrintLevel = 1;
    Tempest::args.bPrintCandidates = 0;
    Tempest::args.deviceOverride = 0;
    Tempest::args.hostMemOverride = 0;
    Tempest::args.deviceMemOverride = 0;
	
    //params
    Tempest::params.sSpectra = 0;
    Tempest::params.sFasta = 0;
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
    Tempest::params.AIons = 0.0;
    Tempest::params.BIons = 0.0;
    Tempest::params.CIons = 0.0;
    Tempest::params.XIons = 0.0;
    Tempest::params.YIons = 0.0;
    Tempest::params.ZIons = 0.0;
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
        // amino acid monoisotopic masses (see constants.h)
        Tempest::params.dMassAA['A'+offset] = A_AA_MASS;
        Tempest::params.dMassAA['B'+offset] = B_AA_MASS;
        Tempest::params.dMassAA['C'+offset] = C_AA_MASS;
        Tempest::params.dMassAA['D'+offset] = D_AA_MASS;
        Tempest::params.dMassAA['E'+offset] = E_AA_MASS;
        Tempest::params.dMassAA['F'+offset] = F_AA_MASS;
        Tempest::params.dMassAA['G'+offset] = G_AA_MASS;
        Tempest::params.dMassAA['H'+offset] = H_AA_MASS;
        Tempest::params.dMassAA['I'+offset] = I_AA_MASS;
        Tempest::params.dMassAA['J'+offset] = J_AA_MASS;
        Tempest::params.dMassAA['K'+offset] = K_AA_MASS;
        Tempest::params.dMassAA['L'+offset] = L_AA_MASS;
        Tempest::params.dMassAA['M'+offset] = M_AA_MASS;
        Tempest::params.dMassAA['N'+offset] = N_AA_MASS;
        Tempest::params.dMassAA['O'+offset] = O_AA_MASS;
        Tempest::params.dMassAA['P'+offset] = P_AA_MASS;
        Tempest::params.dMassAA['Q'+offset] = Q_AA_MASS;
        Tempest::params.dMassAA['R'+offset] = R_AA_MASS;
        Tempest::params.dMassAA['S'+offset] = S_AA_MASS;
        Tempest::params.dMassAA['T'+offset] = T_AA_MASS;
        Tempest::params.dMassAA['U'+offset] = U_AA_MASS;
        Tempest::params.dMassAA['V'+offset] = V_AA_MASS;
        Tempest::params.dMassAA['W'+offset] = W_AA_MASS;
        Tempest::params.dMassAA['X'+offset] = X_AA_MASS;
        Tempest::params.dMassAA['Y'+offset] = Y_AA_MASS;
        Tempest::params.dMassAA['Z'+offset] = Z_AA_MASS;

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
        Tempest::params.unModAA['X'+offset] = 'X';
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
