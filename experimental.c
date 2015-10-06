/*
 *  tempest.c
 *
 *  Created by Brendan Faherty on 12/8/08.
 *  Copyright 2008 Dartmouth College. All rights reserved.
 *
 */

#include "tempest.h"

std::vector<int>   host_iPeakBins;
std::vector<float> host_fPeakInts;

int parse_scan(MSToolkit::Spectrum, double*, int*, int*, int*, float*);
void setup_spectrum(int, int, double, int);
void prepare_peaks(int, float*, int, int, float, long*);
void setup_mod_index(int);

extern void collect_msms_spectra()
{
	char sExpFile[STRING_SIZE];
	FILE *fp;
	long lPeaksBuffer;
	int iSpectrum;
	float *fSpectrum;
	int iParsedPeaks;
	double dPrecursorMass;
	int iPrecursorCharge;
	int iMinBin;
	int iMaxBin;
	float fMaxIntensity;
	int iNumMods;
	char c;
	float fProgress;
	int iMaxMaxBin;
    int iScanNum;

    MSToolkit::MSReader r;
    MSToolkit::Spectrum s;
    int j;

    r.setFilter(MSToolkit::MS2);
    
	// open dta list
	if (r.readFile(args.sSpectra, s) != true) {
		fprintf(stderr, "\nERROR\tUnable to open spectra data file %s: %s\n", args.sSpectra, strerror(errno));
		tempest_exit(EXIT_FAILURE);
	}
	
	// count dtas
	tempest.iNumSpectra = 0;
    iScanNum = 0;
	while ((iScanNum = s.getScanNumber()) > 0) {
        if (iScanNum >= config.minScan && iScanNum <= config.maxScan)
            tempest.iNumSpectra++;
        r.readFile(NULL, s);
	}

	r.readFile(args.sSpectra, s);

	if (args.iPrintLevel) {
		printf(" » Preparing %d MS/MS spectra...     ", tempest.iNumSpectra);
	}

	// Allocate host memory for peak data
	if (0 == (host_iPeakCounts = (int*) calloc(tempest.iNumSpectra,sizeof(int)))) {
		fprintf(stderr, "\nERROR\tUnable allocate host memory for peak counts: %s\n", strerror(errno));
		tempest_exit(EXIT_FAILURE);
	}

	if (0 == (host_lPeakIndices = (long*) malloc (tempest.iNumSpectra*sizeof(long)))) {
		fprintf(stderr, "\nERROR\tUnable allocate host memory for peak indices: %s\n", strerror(errno));
		tempest_exit(EXIT_FAILURE);
	}
	
	// Allocate pinned host memory for candidate data buffers (2 per spectrum)
    //MEA: pin?
    host_cCandidates = (cObj*) malloc(tempest.iNumSpectra * config.iCandidateBufferSize * sizeof(cObj));

	// initialize
	iSpectrum = 0;
    iScanNum = 0;
	tempest.lNumMS2Peaks = 0;
	tempest.fMinPrecursorMass = FLT_MAX;
	tempest.fMaxPrecursorMass = 0.0f;
	iMaxMaxBin = 0;
	
	// Process scans
	while ((iScanNum = s.getScanNumber()) > 0) {
        if (iScanNum < config.minScan || iScanNum > config.maxScan) {
            r.readFile(NULL, s);
            continue;
        }
		// get scan info
		iParsedPeaks = parse_scan(s, &dPrecursorMass, &iPrecursorCharge, &iMinBin, &iMaxBin, &fMaxIntensity);
		if (iParsedPeaks == 0) {
			tempest.iNumSpectra--;
			fprintf(stderr, "WARNING! No peaks loaded from scan #%d.\n", iScanNum);
            r.readFile(NULL, s);
			continue;
		}
        else if (iParsedPeaks == -1) {
            tempest.iNumSpectra--;
            //MEA: warn
			//fprintf(stderr, "WARNING! Precursor mass outside of digest mass range in scan #%d.\n", iScanNum);
            r.readFile(NULL, s);
			continue;
        }

        host_iPeakCounts[iSpectrum] = iParsedPeaks;
        host_lPeakIndices[iSpectrum] = tempest.lNumMS2Peaks;
        tempest.lNumMS2Peaks += iParsedPeaks;
        
		// setup the MS/MS object
		setup_spectrum(iSpectrum, iScanNum, dPrecursorMass, iPrecursorCharge);

		// update min/max precursor masses and max max bin
	 	if ((float) dPrecursorMass < tempest.fMinPrecursorMass) {
	     	tempest.fMinPrecursorMass = (float) dPrecursorMass;
	 	}
		
	 	if ((float) dPrecursorMass > tempest.fMaxPrecursorMass) {
	     	tempest.fMaxPrecursorMass = (float) dPrecursorMass;
	 	}

	 	if (iMaxBin > iMaxMaxBin) {
	 		iMaxMaxBin = iMaxBin;
	 	}

	 	// increment
		iSpectrum++;

		if (args.iPrintLevel) {
			fProgress = (100.0f * iSpectrum) / tempest.iNumSpectra;
			if (fProgress - floor(fProgress) < 0.01f) printf("\b\b\b\b%*d%%", 3, (int) fProgress);
			fflush(0);
		}
        r.readFile(NULL, s);
	}

	// close reader
    //nothing to do?

	if (args.iPrintLevel) {
		printf("\n");
	}
	
	// Setup modified boolean lists
	for (iNumMods = 1; iNumMods <= params.iModificationsMax; iNumMods++) {
		setup_mod_index(iNumMods);
	}

    if (params.bCrossCorrelation)
        tempest.iNumMS2Bins = iMaxMaxBin + tempest.iCrossCorrelationWidth + 1;
    else
        tempest.iNumMS2Bins = iMaxMaxBin + 1;
}

int parse_scan(MSToolkit::Spectrum spectrum, double *dPrecursorMass, int* iPrecursorCharge, int* iMinBin, int* iMaxBin, float* fMaxIntensity) {

	static int iBin;
	static int iNumPeaks;
	static int iCharge;
	static int iRange;
	static float fMz;
	static float fIntensity;
	static FILE* fp;
	static bool bSkip;
    static int j;

    if (spectrum.getCharge() < 0) {
        fprintf(stderr, "\nERROR\tUnable to parse precursor charge from scan #%d.\n", spectrum.getScanNumber());
		tempest_exit(EXIT_FAILURE);
    }
    
    *iPrecursorCharge = spectrum.getCharge();
    
    if (spectrum.getMZ() < 0) {
        fprintf(stderr, "\nERROR\tUnable to parse precursor M/Z from scan #%d.\n", spectrum.getScanNumber());
		tempest_exit(EXIT_FAILURE);
    }
    
    //Mass is M+H
    *dPrecursorMass = spectrum.getMZ() * *iPrecursorCharge - (*iPrecursorCharge-1)*PROTON_MASS;
	
	// validate mass and charge from top line
	if (std::isnan(*dPrecursorMass) || *dPrecursorMass < 1) {
		fprintf(stderr, "\nERROR\tInvalid precursor mass (%f) in scan #%d.\n", *dPrecursorMass, spectrum.getScanNumber());
		return 0;
	}
	
	if (std::isnan(*iPrecursorCharge) || *iPrecursorCharge < 1) {
		fprintf(stderr, "\nERROR\tInvalid precursor charge (%d) in scan #%d.\n", *iPrecursorCharge, spectrum.getScanNumber());
		return 0;
	}

    // return -1 if precursorMass is outside mass range
    if (*dPrecursorMass < MIN_DIGEST_MASS || *dPrecursorMass > MAX_DIGEST_MASS) {
        return -1;
    }
        
	// initialize
	iNumPeaks = spectrum.size();
	*iMinBin = INT_MAX;
	*iMaxBin = 0;
	*fMaxIntensity = 0;

    std::map<int,float> bin2intensity;
    
	// parse peaks
	for (j=0; j<iNumPeaks; j++) {
        fMz = spectrum.at(j).mz;
        fIntensity = spectrum.at(j).intensity;
		// check precursors
		if (params.bRemovePrecursors) {
			for (iCharge=1; iCharge <= *iPrecursorCharge; iCharge++) {
				bSkip = (fabs(fMz - (*dPrecursorMass + iCharge*H_MASS) / iCharge) < params.fFragmentTolerance);
				if (bSkip) break;
			}
			
			if (bSkip) {
				bSkip = 0;
				continue;
			}
		}

		// check clear mz ranges
		for (iRange=0; iRange<params.iNumRemoveMzRanges; iRange++) {
			bSkip = (fMz >= params.fRemoveMzRanges[2*iRange]   - params.fFragmentTolerance &&
				     fMz <= params.fRemoveMzRanges[2*iRange+1] + params.fFragmentTolerance);
				
			if (bSkip) break;
		}
		
		if (bSkip) {
			bSkip = 0;
			continue;
		}

		// get bin
		iBin = (int) (fMz / params.fFragmentTolerance + 0.5f);
		
		// check bin
		if (iBin < 0) {
			fprintf(stderr, "\ninvalid peak (mz: %f, bin: %d)\n", fMz, iBin);
			tempest_exit(EXIT_FAILURE);
		}

        if (fIntensity > bin2intensity[iBin])
            bin2intensity[iBin] = fIntensity;
        
		// update max bin and max intensity
		if (iBin < *iMinBin) *iMinBin = iBin;
		if (iBin > *iMaxBin) *iMaxBin = iBin;
		if (fIntensity > *fMaxIntensity) *fMaxIntensity = fIntensity;
	}

    //setup regioning
	float binsPerRegion = (float)(*iMaxBin-*iMinBin) / MSMS_NORM_REGIONS;
	if (binsPerRegion == 0) binsPerRegion = 1;
	
	//setup noise removal (squared because we sqrt *after* thresholding)
	float fIntensityCutoff = *fMaxIntensity * MSMS_NOISE_CUTOFF * MSMS_NOISE_CUTOFF;
	
	// get max intensity in each region
    float fRegionMaxes[MSMS_NORM_REGIONS];
    for (int iRegion=0; iRegion<MSMS_NORM_REGIONS; iRegion++) {
		fRegionMaxes[iRegion] = 0.0;
	}
    for (std::map<int,float>::iterator it=bin2intensity.begin(); it!=bin2intensity.end(); ++it) {
        iBin = it->first;
        fIntensity = it->second;
        int iRegion = (int)((iBin-*iMinBin)/MSMS_NORM_REGIONS / binsPerRegion * MSMS_NORM_REGIONS);
        if (iRegion == MSMS_NORM_REGIONS)
            iRegion--;
        if (fIntensity > fRegionMaxes[iRegion])
			fRegionMaxes[iRegion] = fIntensity;
    }

    int parsedPeaks = 0;
    for (std::map<int,float>::iterator it=bin2intensity.begin(); it!=bin2intensity.end(); ++it) {
        iBin = it->first;
        fIntensity = it->second;
        if (fIntensity < fIntensityCutoff)
            continue;
        int iRegion = (int)((iBin-*iMinBin)/MSMS_NORM_REGIONS / binsPerRegion * MSMS_NORM_REGIONS);
        if (iRegion == MSMS_NORM_REGIONS)
            iRegion--;
        host_iPeakBins.push_back(iBin);
        host_fPeakInts.push_back(sqrt(fIntensity / fRegionMaxes[iRegion]) * (float)PRIMARY_INTENSITY);
        parsedPeaks++;
    }
    
	return parsedPeaks;
}

int parse_mzxml() {
	return 0;
}

/*
 * Allocate and initialize an MS/MS scan eObj.
 *
 * 		sName 				dta filename, scan number, or other identifier
 * 		dPrecursorMass		precursor mass
 * 		iPrecursorCharge 	precursor charge
 * 		lIndex 				global index of the MS/MS spectrum (in global data arrays)
 */

void setup_spectrum(int iSpectrumIndex, int iScanNum, double dPrecursorMass, int iPrecursorCharge) {
	static eObj *e;
	static eObj *p;
	static int iMassIndex;
    
	// allocate memory
	if (0 == (e = (eObj*) malloc(sizeof(eObj)))) {
		fprintf(stderr, "\nERROR\tUnable to allocate memory for scan object\n");
		tempest_exit(EXIT_FAILURE);
	}

	// initialize
    char tmp[256];
    sprintf(tmp, "%d", iScanNum);
    e->sName = strdup(tmp);
	e->lIndex = iSpectrumIndex;
	e->dPrecursorMass = dPrecursorMass;
	e->iPrecursorCharge = iPrecursorCharge;
	
	e->iNumBufferedCandidates = 0;
	e->iNumCandidates = 0;

	e->pCandidateBufferFill = &host_cCandidates[config.iCandidateBufferSize * e->lIndex];

    e->clEventSent = clCreateUserEvent(clContext, NULL);
    clSetUserEventStatus(e->clEventSent, 0);

	e->next = 0;

	// get mass index
	iMassIndex = (int) floor((float) e->dPrecursorMass / tempest.fPrecursorBinWidth);
		
	// store scan in lookup table
	if (0 == eScanIndex[iMassIndex]) {
		eScanIndex[iMassIndex] = e;
	}
	else {
		p = eScanIndex[iMassIndex];
		while (p->next) p = p->next;
		p->next = e;
	}

	// set unmodified mass booleans
    bScanIndex[0][iMassIndex] = 1;
    if (iMassIndex > 0)
        bScanIndex[0][iMassIndex-1] = 1;
    if (iMassIndex < tempest.iNumPrecursorBins-1)
        bScanIndex[0][iMassIndex+1] = 1;
}

/*
 * Set the MS/MS scan mass index for the given number of mods from the "previous" index.
 */
//MEA: what is this good for?
void setup_mod_index(int iNumMods) {
	int iBin;
	int iModBin;
	int iMod;

	for (iBin=0; iBin<tempest.iNumPrecursorBins; iBin++) {
		
		//skip empty bins
		if (0 == bScanIndex[iNumMods-1][iBin]) continue;
		
		// aa mods
		for (iMod=0; iMod<params.iNumMods; iMod++) {
			if (params.tMods[iMod].cSymbol) {
				iModBin = (int) (iBin - params.tMods[iMod].dMassDiff / tempest.fPrecursorBinWidth);
				bScanIndex[iNumMods][iModBin] = 1;
			}
		}

		// nterm
		if (params.cVariableNtermSymbol) {
			iModBin = (int) (iBin - params.dVariableNtermMassDiff / tempest.fPrecursorBinWidth);
			bScanIndex[iNumMods][iModBin] = 1;
		}

		// cterm
		if (params.cVariableCtermSymbol) {
			iModBin = (int) (iBin - params.dVariableCtermMassDiff / tempest.fPrecursorBinWidth);
			bScanIndex[iNumMods][iModBin] = 1;
		}
	}
}

