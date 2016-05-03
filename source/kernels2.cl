#define STRING_SIZE 256
#define CASE_SHIFT 32
#define MIN_SCORE -1000.0f

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
    cObj   cCandidate;
    int    iNumOccurrences;
    float  fScore;
} mObj;

//Mass delta and score weighting for a neutral loss event
typedef struct nlValue {
    int   numSites;
    char  hasAA[256]; //boolean
    char  modNum; //values from 0 to 5
    float massDelta;
    float weighting;
} nlValue;

//=================================================================================================
// Kernel Functions
//=================================================================================================

__kernel void cl_build(__global float* fSpectrum, int iPeakCount, __global int* iPeakBins, __global float* fPeakInts, long lPeakOffset)
{
    int iBin;
    float intensity;
    
    int iThread = get_global_id(0);
    
    iPeakBins += lPeakOffset;
    fPeakInts += lPeakOffset;

    if (iThread < iPeakCount) {
        iBin = iPeakBins[iThread];
        intensity = fPeakInts[iThread];
        fSpectrum[iBin] = intensity;
    }    
}

__kernel void cl_transform(__global float* fSpectrum, int NUM_BINS)
{
    int iStart = get_global_id(0);
    if (iStart >= NUM_BINS) return;

    float fIntensity = 0.0f;

    // Calculate the initial fIntensity for each thread
    int windowStart = max(iStart-XCORR_TRANSFORM_WIDTH, 0);
    int windowEnd = min(iStart+XCORR_TRANSFORM_WIDTH, NUM_BINS-1);
    for (int iBin = windowStart; iBin <= windowEnd; iBin++) {
        fIntensity += fSpectrum[iBin];
    }

    barrier(CLK_GLOBAL_MEM_FENCE);
    
    fSpectrum[iStart] -= (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
}

__kernel void cl_score(int iPrecursorCharge, int iNumCandidates, __global cObj* gpu_cCandidates, __global float* gpu_fScores, __global float* gpu_fSpectra, long lSpectrumOffset, int NUM_BINS, __constant float* MASS_AA, __constant float* NTERM_MOD_MASSES, constant float* CTERM_MOD_MASSES, __constant nlValue* nlValuesNterm, __constant nlValue* nlValuesCterm, __constant nlValue* nlValuesAA)
{

    if (get_global_id(0) >= iNumCandidates) {
        gpu_fScores[get_global_id(0)] = 0.0f;
        return;
    }
  
    __global float* fSpectrum = gpu_fSpectra + lSpectrumOffset;
    cObj cCandidate = gpu_cCandidates[get_global_id(0)];
    
    int iBin;
    float nTermFragMass;
    float baseMass = cCandidate.fPeptideMass;
    float fScore = 0.0f;

#if NUM_AA_NL
    //if !NUM_AA_NL, these won't be declared and won't take up register space.
    int nNLTotal = 0;
    int cNLTotal = 0;
    int nNLCount[NUM_AA_NL];
    int cNLCount[NUM_AA_NL];
#endif
    
#if NUM_AA_NL
    //zero count arrays
    for (int i=0; i<NUM_AA_NL; i++) {
        nNLCount[i] = 0;
        cNLCount[i] = 0;
        //initialize aa counts to y fragment
        for(int iFragment=0; iFragment<cCandidate.iPeptideLength; iFragment++) {
            if (nlValuesAA[i].hasAA[cCandidate.sPeptide[iFragment]]) {
                cNLCount[i] += 1;
                cNLTotal += 1;
            }
        }
    }
#endif
    
    // Base Masses
    nTermFragMass = MASS_PEPTIDE_NTERM;
    if (cCandidate.cBefore == '-') nTermFragMass += MASS_PROTEIN_NTERM;
    if (cCandidate.ntermMod) nTermFragMass += NTERM_MOD_MASSES[cCandidate.ntermMod-1];

    for (int iFragment=0; iFragment<cCandidate.iPeptideLength-1; iFragment++) {
        unsigned char aa = cCandidate.sPeptide[iFragment];
        //update phos mods
#if NUM_AA_NL
        if (cNLTotal) {
            for (int i=0; i<NUM_AA_NL; i++) {
                if (nlValuesAA[i].hasAA[aa]) {
                    nNLCount[i] += 1;
                    cNLCount[i] -= 1;
                    nNLTotal += 1;
                    cNLTotal -= 1;
                }
            }
        }
#endif

        //calculate the fragment mass
        nTermFragMass += MASS_AA[aa];
        //baseMass -= MASS_AA[aa];

        for (int iCharge=1; iCharge<iPrecursorCharge; iCharge++) {
            
            if (A_IONS != 0) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * A_IONS;
            }
            if (B_IONS != 0) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * B_IONS;
            }
            if (C_IONS != 0) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * C_IONS;
            }

            if (X_IONS != 0) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * X_IONS;
            }
            if (Y_IONS != 0) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1)) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * Y_IONS;
            }
            if (Z_IONS != 0) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < NUM_BINS)
                    fScore += fSpectrum[iBin] * Z_IONS;
            }
            
#if NUM_AA_NL
            if (cNLTotal || nNLTotal){
                for (int i=0; i<NUM_AA_NL; i++) {
                    if (nNLCount[i]) {
                        if (A_IONS != 0) {
                            iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * A_IONS * nlValuesAA[i].weighting;
                        }
                        if (B_IONS != 0) {
                            iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * B_IONS * nlValuesAA[i].weighting;
                        }
                        if (C_IONS != 0) {
                            iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * C_IONS * nlValuesAA[i].weighting;
                        }
                    }
                    if (cNLCount[i]) {
                        if (X_IONS != 0) {
                            iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * X_IONS * nlValuesAA[i].weighting;
                        }
                        if (Y_IONS != 0) {
                            iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * Y_IONS * nlValuesAA[i].weighting;
                        }
                        if (Z_IONS != 0) {
                            iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS - nlValuesAA[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                            if (0 <= iBin && iBin < NUM_BINS)
                                fScore += fSpectrum[iBin] * Z_IONS * nlValuesAA[i].weighting;
                        }
                    }
                }
            }
#endif                   

#if NUM_NTERM_NL
            for (int i=0; i<NUM_NTERM_NL; i++) {
                if (nlValuesNterm[i].modNum != cCandidate.ntermMod)
                    continue;
                if (A_IONS != 0) {
                    iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS - nlValuesNterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * A_IONS * nlValuesNterm[i].weighting;
                }
                if (B_IONS != 0) {
                    iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - nlValuesNterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * B_IONS * nlValuesNterm[i].weighting;
                }
                if (C_IONS != 0) {
                    iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS - nlValuesNterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * C_IONS * nlValuesNterm[i].weighting;
                }
            }
#endif
#if NUM_CTERM_NL
            for (int i=0; i<NUM_CTERM_NL; i++) {
                if (nlValuesCterm[i].modNum != cCandidate.ctermMod)
                    continue;
                if (X_IONS != 0) {
                    iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS - nlValuesCterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * X_IONS * nlValuesCterm[i].weighting;
                }
                if (Y_IONS != 0) {
                    iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - nlValuesCterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * Y_IONS * nlValuesCterm[i].weighting;
                }
                if (Z_IONS != 0) {
                    iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS - nlValuesCterm[i].massDelta) / iCharge / FRAGMENT_TOLERANCE + FRAGMENT_BIN_OFFSET);
                    if (0 <= iBin && iBin < NUM_BINS)
                        fScore += fSpectrum[iBin] * Z_IONS * nlValuesCterm[i].weighting;
                }
            }
#endif
        }
    }                               

    gpu_fScores[get_global_id(0)] = fScore*PRIMARY_INTENSITY/10000.0f;
    
}

__kernel void cl_reduce_scores_sequential(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset)
{
    int iThread = get_global_id(0);
    if (iThread > 0) return;

    mPSMs += lPSMsOffset;

    //__local int topNMatches[NUM_INTERNAL_PSMS];
    float topNScores[NUM_INTERNAL_PSMS];
    float score;
    float minScore;
    minScore = mPSMs[NUM_INTERNAL_PSMS-1].fScore;

    for (int i=0; i<NUM_INTERNAL_PSMS; i++) {
        topNScores[i] = mPSMs[i].fScore;
    }
    
    for (int i=0; i<iNumScores; i++) {
        score = fCandidateScores[i];
        if (score < minScore)
            continue;
        
        for (int j=0; j<NUM_INTERNAL_PSMS; j++) {
            if (score > topNScores[j]) {
                //shift right and insert
                for (int k=NUM_INTERNAL_PSMS-1; k>j; k--) {
                    topNScores[k] = topNScores[k-1];
                    mPSMs[k] = mPSMs[k-1];
                }
                topNScores[j] = score;
                minScore = topNScores[NUM_INTERNAL_PSMS-1];

                mPSMs[j].cCandidate = cCandidates[i];
                mPSMs[j].fScore = score;
                mPSMs[j].iNumOccurrences = !cCandidates[i].decoy; //don't count decoys as occurrences
                break;
            }
            else if (score == topNScores[j]) {
                int aa=0;
                while (aa<mPSMs[j].cCandidate.iPeptideLength && (mPSMs[j].cCandidate.sPeptide[aa] == cCandidates[i].sPeptide[aa]))
                    aa++;
                if (mPSMs[j].cCandidate.sPeptide[aa] == cCandidates[i].sPeptide[aa]) {
                    //duplicate found: report lowest non-decoy protein ID and increment occurrences
                    if (!(!mPSMs[j].cCandidate.decoy && cCandidates[i].decoy)) {
                        if (cCandidates[i].iProtein < mPSMs[j].cCandidate.iProtein) { 
                            mPSMs[j].cCandidate.iProtein = cCandidates[i].iProtein;
                            mPSMs[j].cCandidate.cBefore = cCandidates[i].cBefore;
                            mPSMs[j].cCandidate.cAfter = cCandidates[i].cAfter;
                        }
                        mPSMs[j].iNumOccurrences += !cCandidates[i].decoy;
                    }
                    break;
                }
            }
        }
    }
}

__kernel void cl_reduce_scores_parallel(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset, __local int* iMatches, __local float* fScores)
{
    int iThread = get_local_id(0);
    __local int psmInd;
    __local cObj topCandidate;
    psmInd = 0;
    
    mPSMs += lPSMsOffset;

    // thread, scores per thread, and thread index
    int iScoresPerThread = iNumScores/get_local_size(0);
    int iThreadIndex = iThread*iScoresPerThread;
    
    // Part One: Each thread reduces a portion of the scores
    while (psmInd < NUM_INTERNAL_PSMS) {
        // initialize local memory
        iMatches[iThread] = iThreadIndex;
        fScores[iThread] = fCandidateScores[iThreadIndex];

        // loop through scores for this thread
        for (int i = iThreadIndex+1; i < (iThread + 1) * iScoresPerThread; i++) {
            float f = fCandidateScores[i];
            if (f > fScores[iThread]) {
                iMatches[iThread] = i;
                fScores[iThread] = f;
            }                               
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (int offset=get_local_size(0)/2; offset>0; offset/=2) {
            if (iThread<offset) {
                //gpu_reduce_score(fScores, iMatches, iThread, i);
                if (fScores[iThread] < fScores[iThread + offset]) {
                    iMatches[iThread] = iMatches[iThread + offset];
                    fScores[iThread] = fScores[iThread + offset];
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // thread 0 holds the top PSM
        if (iThread != 0)
            goto sync;
        
        if (fScores[0] < mPSMs[NUM_INTERNAL_PSMS-1].fScore || fScores[0] <= 0) {
            psmInd = NUM_INTERNAL_PSMS; //break from loops
            goto sync;
        }

        topCandidate = cCandidates[iMatches[0]];
        //new best score
        while (psmInd<NUM_INTERNAL_PSMS && fScores[0]<mPSMs[psmInd].fScore)
            psmInd++;
        for (int i=psmInd; i<NUM_INTERNAL_PSMS; i++) {
            if (fScores[0] > mPSMs[i].fScore) {
                //shift right
                for (int j=NUM_INTERNAL_PSMS-1; j>i; j--)
                    mPSMs[j] = mPSMs[j-1];

                // update ith PSM
                mPSMs[i].cCandidate = topCandidate;
                mPSMs[i].fScore = fScores[0];
                mPSMs[i].iNumOccurrences = !topCandidate.decoy; //don't count decoys as occurrences

                //currentScore = fScores[0];
                break;
            }
            else if (fScores[0] == mPSMs[i].fScore) {
                int aa=0;
                while (aa<mPSMs[i].cCandidate.iPeptideLength && (mPSMs[i].cCandidate.sPeptide[aa] == topCandidate.sPeptide[aa]))
                    aa++;
                if (mPSMs[i].cCandidate.sPeptide[aa] == topCandidate.sPeptide[aa]) {
                    //duplicate found: report lowest non-decoy protein ID and increment occurrences
                    if (!(!mPSMs[i].cCandidate.decoy && topCandidate.decoy)) {
                        if (topCandidate.iProtein < mPSMs[i].cCandidate.iProtein) {
                            mPSMs[i].cCandidate.iProtein = topCandidate.iProtein;
                            mPSMs[i].cCandidate.cBefore = topCandidate.cBefore;
                            mPSMs[i].cCandidate.cAfter = topCandidate.cAfter;
                        }
                        mPSMs[i].iNumOccurrences += !topCandidate.decoy;
                    }
                    break;
                }
            }
        }
        fCandidateScores[iMatches[0]] = MIN_SCORE;
        
    sync:
        barrier(CLK_LOCAL_MEM_FENCE);
    }
}
