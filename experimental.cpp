#include "tempest.hpp"

int parse_scan(MSToolkit::Spectrum, double*, int*, int*, int*, float*);
void setup_spectrum(int, int, double, int);
void prepare_peaks(int, float*, int, int, float, long*);
//void setup_mod_index(int);

extern void Tempest::collect_msms_spectra()
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
	if (r.readFile(Tempest::args.sSpectra, s) != true) {
		fprintf(stderr, "\nERROR\tUnable to open spectra data file %s: %s\n", Tempest::args.sSpectra, strerror(errno));
		Tempest::tempest_exit(EXIT_FAILURE);
	}
	
	// count dtas
	Tempest::data.iNumSpectra = 0;
    iScanNum = 0;
	while ((iScanNum = s.getScanNumber()) > 0) {
        if (iScanNum >= config.minScan && iScanNum <= config.maxScan)
            Tempest::data.iNumSpectra++;
        r.readFile(NULL, s);
	}

	r.readFile(Tempest::args.sSpectra, s);

	if (Tempest::args.iPrintLevel) {
		printf(" » Preparing %d MS/MS spectra...     ", Tempest::data.iNumSpectra);
	}

	// Allocate host memory for peak data
	if (0 == (Tempest::data.host_iPeakCounts = (int*) calloc(Tempest::data.iNumSpectra,sizeof(int)))) {
		fprintf(stderr, "\nERROR\tUnable allocate host memory for peak counts: %s\n", strerror(errno));
		Tempest::tempest_exit(EXIT_FAILURE);
	}

	if (0 == (Tempest::data.host_lPeakIndices = (long*) malloc (Tempest::data.iNumSpectra*sizeof(long)))) {
		fprintf(stderr, "\nERROR\tUnable allocate host memory for peak indices: %s\n", strerror(errno));
		Tempest::tempest_exit(EXIT_FAILURE);
	}
	
	// Allocate pinned host memory for candidate data buffers (2 per spectrum)
    //MEA: pin?
    //Tempest::data.host_cCandidates = (cObj*) malloc(Tempest::data.iNumSpectra * config.iCandidateBufferSize * sizeof(cObj));

	// initialize
	iSpectrum = 0;
    iScanNum = 0;
	Tempest::data.lNumMS2Peaks = 0;
	Tempest::data.fMinPrecursorMass = FLT_MAX;
	Tempest::data.fMaxPrecursorMass = 0.0f;
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
			Tempest::data.iNumSpectra--;
			fprintf(stderr, "WARNING! No peaks loaded from scan #%d.\n", iScanNum);
            r.readFile(NULL, s);
			continue;
		}
        else if (iParsedPeaks == -1) {
            Tempest::data.iNumSpectra--;
            //MEA: warn
			//fprintf(stderr, "WARNING! Precursor mass outside of digest mass range in scan #%d.\n", iScanNum);
            r.readFile(NULL, s);
			continue;
        }

        Tempest::data.host_iPeakCounts[iSpectrum] = iParsedPeaks;
        Tempest::data.host_lPeakIndices[iSpectrum] = Tempest::data.lNumMS2Peaks;
        Tempest::data.lNumMS2Peaks += iParsedPeaks;
        
		// setup the MS/MS object
		setup_spectrum(iSpectrum, iScanNum, dPrecursorMass, iPrecursorCharge);

		// update min/max precursor masses and max max bin
	 	if ((float) dPrecursorMass < Tempest::data.fMinPrecursorMass) {
	     	Tempest::data.fMinPrecursorMass = (float) dPrecursorMass;
	 	}
		
	 	if ((float) dPrecursorMass > Tempest::data.fMaxPrecursorMass) {
	     	Tempest::data.fMaxPrecursorMass = (float) dPrecursorMass;
	 	}

	 	if (iMaxBin > iMaxMaxBin) {
	 		iMaxMaxBin = iMaxBin;
	 	}

	 	// increment
		iSpectrum++;

		if (Tempest::args.iPrintLevel) {
			fProgress = (100.0f * iSpectrum) / Tempest::data.iNumSpectra;
			if (fProgress - floor(fProgress) < 0.01f) printf("\b\b\b\b%*d%%", 3, (int) fProgress);
			fflush(0);
		}
        r.readFile(NULL, s);
	}

	// close reader
    //nothing to do?

	if (Tempest::args.iPrintLevel) {
		printf("\n");
	}
	
	// Setup modified boolean lists
	// for (iNumMods = 1; iNumMods <= Tempest::params.iModificationsMax; iNumMods++) {
	// 	setup_mod_index(iNumMods);
	// }

    //if (Tempest::params.bCrossCorrelation)
    //    Tempest::data.iNumMS2Bins = iMaxMaxBin + Tempest::data.iCrossCorrelationWidth + 1;
    //else
    Tempest::data.iNumMS2Bins = iMaxMaxBin + 1;
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
		Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    *iPrecursorCharge = spectrum.getCharge();
    
    if (spectrum.getMZ() < 0) {
        fprintf(stderr, "\nERROR\tUnable to parse precursor M/Z from scan #%d.\n", spectrum.getScanNumber());
		Tempest::tempest_exit(EXIT_FAILURE);
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
    if (*dPrecursorMass < Tempest::params.minPeptideMass || *dPrecursorMass > Tempest::params.maxPeptideMass) {
        return -1;
    }
        
	// initialize
	iNumPeaks = spectrum.size();
	*iMinBin = INT_MAX;
	*iMaxBin = 0;
	*fMaxIntensity = 0;

    std::map<int,float> bin2intensity;

    if (Tempest::params.intensityThreshold > 0)
        for (j=0; j<iNumPeaks; j++)
            if (spectrum.at(j).intensity > *fMaxIntensity)
                *fMaxIntensity = spectrum.at(j).intensity;

    float intensityCutoff = *fMaxIntensity * Tempest::params.intensityThreshold;
    
	// parse peaks
	for (j=0; j<iNumPeaks; j++) {
        fMz = spectrum.at(j).mz;
        fIntensity = spectrum.at(j).intensity;

        if (fIntensity < intensityCutoff)
            continue;
        
		// check precursors
		if (Tempest::params.bRemovePrecursors) {
			for (iCharge=1; iCharge <= *iPrecursorCharge; iCharge++) {
				bSkip = (fabs(fMz - (*dPrecursorMass + iCharge*PROTON_MASS) / iCharge) < Tempest::params.fFragmentTolerance);
				if (bSkip) break;
			}
			
			if (bSkip) {
				bSkip = 0;
				continue;
			}
		}

		// check clear mz ranges
		for (iRange=0; iRange<Tempest::params.iNumRemoveMzRanges; iRange++) {
			bSkip = (fMz >= Tempest::params.fRemoveMzRanges[2*iRange]   - Tempest::params.fFragmentTolerance &&
				     fMz <= Tempest::params.fRemoveMzRanges[2*iRange+1] + Tempest::params.fFragmentTolerance);
				
			if (bSkip) break;
		}
		
		if (bSkip) {
			bSkip = 0;
			continue;
		}

		// get bin
		iBin = (int) (fMz / Tempest::params.fFragmentTolerance + Tempest::params.fragmentBinOffset);

        // check bin
		if (iBin < 0) {
			fprintf(stderr, "\ninvalid peak (mz: %f, bin: %d)\n", fMz, iBin);
			Tempest::tempest_exit(EXIT_FAILURE);
		}

        //store intensity in primary bin
        bin2intensity[iBin] += fIntensity;
        
		// update max bin and max intensity
		if (iBin < *iMinBin) *iMinBin = iBin;
		if (iBin > *iMaxBin) *iMaxBin = iBin;
		if (fIntensity > *fMaxIntensity) *fMaxIntensity = fIntensity;

        //store intensity in flanking bins
        if (Tempest::params.flankingIntensity > 0) {
            //left flanking bin
            if (iBin > 0) {
                bin2intensity[iBin-1] += fIntensity * Tempest::params.flankingIntensity;
                if (iBin-1 < *iMinBin) *iMinBin = iBin-1;
            }
            //right flanking bin
            bin2intensity[iBin+1] += fIntensity * Tempest::params.flankingIntensity;
            if (iBin+1 > *iMaxBin) *iMaxBin = iBin+1;
        }
	}
    
    //setup regioning
	float binsPerRegion = (float)(*iMaxBin-*iMinBin) / Tempest::params.numNormRegions;
	if (binsPerRegion == 0) binsPerRegion = 1;
	
	// get max intensity in each region
    float fRegionMaxes[Tempest::params.numNormRegions];
    for (int iRegion=0; iRegion<Tempest::params.numNormRegions; iRegion++) {
		fRegionMaxes[iRegion] = 0.0;
	}
    for (std::map<int,float>::iterator it=bin2intensity.begin(); it!=bin2intensity.end(); ++it) {
        iBin = it->first;
        fIntensity = it->second;
        int iRegion = (int)((iBin-*iMinBin)/Tempest::params.numNormRegions / binsPerRegion * Tempest::params.numNormRegions);
        if (iRegion == Tempest::params.numNormRegions)
            iRegion--;
        if (fIntensity > fRegionMaxes[iRegion])
			fRegionMaxes[iRegion] = fIntensity;
    }

    int parsedPeaks = 0;
    for (std::map<int,float>::iterator it=bin2intensity.begin(); it!=bin2intensity.end(); ++it) {
        iBin = it->first;
        fIntensity = it->second;
        int iRegion = (int)((iBin-*iMinBin)/Tempest::params.numNormRegions / binsPerRegion * Tempest::params.numNormRegions);
        if (iRegion == Tempest::params.numNormRegions)
            iRegion--;
        Tempest::data.host_iPeakBins.push_back(iBin);
        //Tempest::data.host_fPeakInts.push_back(pow(fIntensity / fRegionMaxes[iRegion], (float)1/2) * (float)PRIMARY_INTENSITY);
        Tempest::data.host_fPeakInts.push_back(sqrt(fIntensity / fRegionMaxes[iRegion]) * (float)PRIMARY_INTENSITY);
        parsedPeaks++;
    }
    
    /*
    //setup regioning
	float binsPerRegion = (float)Tempest::params.numNormRegions/Tempest::params.fFragmentTolerance;
    if (binsPerRegion == 0) binsPerRegion = 1;
    float halfRegionWidth = binsPerRegion/2;
	
    int parsedPeaks = 0;
    for (std::map<int,float>::iterator it=bin2intensity.begin(); it!=bin2intensity.end(); ++it) {
        iBin = it->first;
        fIntensity = it->second;
        float regionMax = 0;
        for (std::map<int,float>::iterator it2=bin2intensity.begin(); it2!=bin2intensity.end(); ++it2) {
            int iBin2 = it2->first;
            float fIntensity2 = it2->second;
            if (fIntensity2 > regionMax && abs(iBin-iBin2) < halfRegionWidth)
                regionMax = fIntensity2;
            if (iBin2-iBin > halfRegionWidth)
                break;
        }       
        Tempest::data.host_iPeakBins.push_back(iBin);
        Tempest::data.host_fPeakInts.push_back(sqrt(sqrt(fIntensity / regionMax)) * (float)PRIMARY_INTENSITY);
        parsedPeaks++;
    }
    */
    
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
		Tempest::tempest_exit(EXIT_FAILURE);
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

	//e->candidateBuffer = &Tempest::data.host_cCandidates[config.iCandidateBufferSize * e->lIndex];

	e->next = 0;

    // store in 1D vector of scans
    Tempest::data.eScans.push_back(e);
    
	// get mass index
	iMassIndex = (int) floor((float) e->dPrecursorMass / Tempest::data.fPrecursorBinWidth);
		
	// store scan in lookup table
	if (0 == Tempest::data.eScanIndex[iMassIndex]) {
		Tempest::data.eScanIndex[iMassIndex] = e;
	}
	else {
		p = Tempest::data.eScanIndex[iMassIndex];
		while (p->next) p = p->next;
		p->next = e;
	}

	// // set unmodified mass booleans
    // bScanIndex[0][iMassIndex] = 1;
    // if (iMassIndex > 0)
    //     bScanIndex[0][iMassIndex-1] = 1;
    // if (iMassIndex < Tempest::data.iNumPrecursorBins-1)
    //     bScanIndex[0][iMassIndex+1] = 1;
}

/*
 * Set the MS/MS scan mass index for the given number of mods from the "previous" index.
 */
//MEA: what is this good for?
/*
void setup_mod_index(int iNumMods) {
	int iBin;
	int iModBin;
	int iMod;

	for (iBin=0; iBin<Tempest::data.iNumPrecursorBins; iBin++) {
		
		//skip empty bins
		if (0 == bScanIndex[iNumMods-1][iBin]) continue;
		
		// aa mods
		for (iMod=0; iMod<Tempest::params.iNumMods; iMod++) {
			if (Tempest::params.tMods[iMod].cSymbol) {
				iModBin = (int) (iBin - Tempest::params.tMods[iMod].dMassDiff / Tempest::data.fPrecursorBinWidth);
				bScanIndex[iNumMods][iModBin] = 1;
			}
		}

		// nterm
		if (Tempest::params.cVariableNtermSymbol) {
			iModBin = (int) (iBin - Tempest::params.dVariableNtermMassDiff / Tempest::data.fPrecursorBinWidth);
			bScanIndex[iNumMods][iModBin] = 1;
		}

		// cterm
		if (Tempest::params.cVariableCtermSymbol) {
			iModBin = (int) (iBin - Tempest::params.dVariableCtermMassDiff / Tempest::data.fPrecursorBinWidth);
			bScanIndex[iNumMods][iModBin] = 1;
		}
	}
}
*/
