#define STRING_SIZE 256
#define CASE_SHIFT 32

//MEA: include tempest.h instead of copying these?
// scoring
#define MIN_SCORE -1000.0f

// Candidate Peptide
typedef struct cobject {
    int    iProtein;
    float  fPeptideMass;

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
    
    char   ntermMod; //values from 0 to 5
    char   ctermMod; //values from 0 to 5
    char   cBefore;
    char   cAfter;

    unsigned char sPeptide[MAX_PEPTIDE_LENGTH];
} mObj;

//Mass delta and score weighting for a neutral loss event
typedef struct nlValue {
    int   numSites;
    char  hasAA[256]; //boolean
    char  modNum; //values from 0 to 5
    float massDelta;
    float weighting;
} nlValue;

void gpu_reduce_score(__local float* fScores, __local int* iMatches, int iThread, int iThreadOffset);
bool gpu_duplicate(__global cObj *cCandidate, float fCandidateScore, __global mObj* mPSM);
bool duplicate(__global cObj* cCandidate1, __global cObj* cCandidate2);
bool gpu_backbone_match(__global cObj *cCandidate, __global mObj* mPSM);
bool backbone_match(__global cObj* cCandidate1, __global cObj* cCandidate2);
bool gpu_strcmp(__global const unsigned char * s1, __global const unsigned char * s2, int peptideLength);
bool gpu_strcmpi(__global const unsigned char * s1, __global const unsigned char * s2, int peptideLength);

//=================================================================================================
// Global Kernel Functions
//=================================================================================================


__kernel void gpu_build(__global float* fSpectrum, int GPU_NUM_BINS, int iPeakCount, __global int* iPeakBins, __global float* fPeakInts, long lPeakOffset)
{
    //printf("gpu_build\n");
    int iBin, iThread;
    float intensity;
    //iThread = blockIdx.x * blockDim.x + threadIdx.x;
    iThread = get_global_id(0);
    
    iPeakBins += lPeakOffset;
    fPeakInts += lPeakOffset;

    if (iThread < iPeakCount) {
        iBin = iPeakBins[iThread];
        intensity = fPeakInts[iThread];
        fSpectrum[iBin] = intensity;
    }

    if (FLANKING_INTENSITY) {
        barrier(CLK_GLOBAL_MEM_FENCE);
        if (iThread < iPeakCount) {
            if (iBin > 0)
                intensity += fSpectrum[iBin-1]*FLANKING_INTENSITY;
            if (iBin < GPU_NUM_BINS-1)
                intensity += fSpectrum[iBin+1]*FLANKING_INTENSITY;
        }
        barrier(CLK_GLOBAL_MEM_FENCE);

        if (iThread < iPeakCount)
            fSpectrum[iBin] = intensity;
    }
    
}

/*
 * Kernel for transforming a single observed spectrum in global memory.
 * 
 *      fScratch        Built (rasterized) spectrum, prior to transformation.
 *      fSpectrum       Memory for transformed spectrum.
 */

// __kernel void gpu_transform_scratch(__global float* fScratch, long lScratchOffset, __global float* fSpectrum, long lSpectrumOffset)
// {

//     int iStart = get_global_id(0);
//     if (iStart >= GPU_NUM_BINS) return;

//     fScratch += lScratchOffset;
//     fSpectrum += lSpectrumOffset;

//     float fIntensity = 0.0f;

//     // Calculate the initial fIntensity for each thread
//     for (int iBin = iStart-XCORR_TRANSFORM_WIDTH; iBin <= iStart+XCORR_TRANSFORM_WIDTH; iBin++) {
//         fIntensity += fScratch[iBin];
//     }
    
//     fSpectrum[iStart] = (float) fScratch[iStart] - (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
        
//     //update the intensity
//     //fIntensity -= fScratch[iStart-XCORR_TRANSFORM_WIDTH];
//     //fIntensity += fScratch[iStart+XCORR_TRANSFORM_WIDTH+1];
// }

__kernel void gpu_transform(__global float* fSpectrum, int GPU_NUM_BINS)
{

    int iStart = get_global_id(0);
    if (iStart >= GPU_NUM_BINS) return;

    float fIntensity = 0.0f;

    // Calculate the initial fIntensity for each thread
    int windowStart = max(iStart-XCORR_TRANSFORM_WIDTH, 0);
    int windowEnd = min(iStart+XCORR_TRANSFORM_WIDTH, GPU_NUM_BINS-1);
    for (int iBin = windowStart; iBin <= windowEnd; iBin++) {
        fIntensity += fSpectrum[iBin];
    }
    
    fSpectrum[iStart] -= (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
        
    //update the intensity
    //fIntensity -= fScratch[iStart-XCORR_TRANSFORM_WIDTH];
    //fIntensity += fScratch[iStart+XCORR_TRANSFORM_WIDTH+1];
}

// __kernel void gpu_transform_local(__global float* fSpectrum, long lSpectrumOffset)
// {

//     int iStart = get_group_id(0);
//     if (iStart >= GPU_NUM_BINS) return;
//     int localIndex = iStart + (get_local_id(0) - XCORR_TRANSFORM_WIDTH/2);
//     if (localIndex >= iStart)
//         localIndex += 1;

//     fSpectrum += lSpectrumOffset;

//     __local float fIntensity;
//     fIntensity = 0.0f;
//     barrier(CLK_LOCAL_MEM_FENCE);

//     if (localIndex >= 0 && localIndex < GPU_NUM_BINS)
//         atomic_xchg(&fIntensity, fIntensity+fSpectrum[localIndex]);

//     barrier(CLK_LOCAL_MEM_FENCE);
//     barrier(CLK_GLOBAL_MEM_FENCE);
          
//     fSpectrum[iStart] -= (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
        
//     //update the intensity
//     //fIntensity -= fScratch[iStart-XCORR_TRANSFORM_WIDTH];
//     //fIntensity += fScratch[iStart+XCORR_TRANSFORM_WIDTH+1];
// }

// __kernel void gpu_transform2(__global float* fScratch, long lScratchOffset, __global float* fSpectrum, long lSpectrumOffset)
// {
//     int iBin, iBinsPerThread, iStart;
//     float fIntensity;

//     //if (blockIdx.x > 0) return;
//     if (get_group_id(0) > 0) return;

//     fScratch += lScratchOffset;
//     fSpectrum += lSpectrumOffset;

//     fIntensity = 0.0f;

//     //iBinsPerThread = GPU_NUM_BINS / blockDim.x;
//     //iStart = threadIdx.x * iBinsPerThread;
//     iBinsPerThread = GPU_NUM_BINS / get_local_size(0);
//     iStart = get_local_id(0) * iBinsPerThread;

//     // Calculate the initial fIntensity for each thread
//     for (iBin = iStart-XCORR_TRANSFORM_WIDTH; iBin <= iStart+XCORR_TRANSFORM_WIDTH; iBin++) {
//         fIntensity += fScratch[iBin];
//     }
    
//     // Loop through the bins for each thread
//     for (iBin=iStart; iBin+1<iStart+iBinsPerThread; iBin++) {
//         //store the new intensity
//         fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
        
//         //update the intensity
//         fIntensity -= fScratch[iBin-XCORR_TRANSFORM_WIDTH];
//         fIntensity += fScratch[iBin+XCORR_TRANSFORM_WIDTH+1];
//     }
    
//     // use the last thread (if necessary) for the rest of the bins
//     //if (threadIdx.x == blockDim.x-1 && fIntensity>0.0f) {
//     if (get_local_id(0) == get_local_size(0)-1 && fIntensity>0.0f) {
//         for (; iBin+1<GPU_NUM_BINS; iBin++) {
//             //store the new intensity
//             fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
            
//             //update the intensity
//             fIntensity -= fScratch[iBin-XCORR_TRANSFORM_WIDTH];
//             fIntensity += fScratch[iBin+XCORR_TRANSFORM_WIDTH+1];
//         }
//     }
    
//     // and get the last one
//     fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(XCORR_TRANSFORM_WIDTH*2+1);
// }

__kernel void gpu_score(int iPrecursorCharge, int iNumCandidates, __global cObj* gpu_cCandidates, __global float* gpu_fScores, __global float* gpu_fSpectra, long lSpectrumOffset, int GPU_NUM_BINS, __constant float* GPU_MASS_AA, __constant float* NTERM_MOD_MASSES, constant float* CTERM_MOD_MASSES, __constant nlValue* nlValuesNterm, __constant nlValue* nlValuesCterm, __constant nlValue* nlValuesAA)
{
    int iThread = get_global_id(0);
    if (iThread >= iNumCandidates) {
        gpu_fScores[iThread] = 0.0f;
        return;
    }

    //if (iThread == 0)
    //    printf("%d\n", get_local_size(0));
  
    __global float* fSpectrum = gpu_fSpectra + lSpectrumOffset;
    //__local cObj cCandidates[BLOCKDIM_SCORE];
    //cCandidates[get_local_id(0)] = gpu_cCandidates[get_group_id(0)*BLOCKDIM_SCORE+get_local_id(0)];

    //__local cObj* cCandidate = &(cCandidates[get_local_id(0)]);
    cObj cCandidate = gpu_cCandidates[iThread];
    
    int iCharge;
    int iBin;
    int iFragment;
    float nTermFragMass;
    float baseMass = cCandidate.fPeptideMass;
    float fScore = 0.0f;
    int i;
    char aa;
    int nNLTotal = 0;
    int cNLTotal = 0;
    int nNLCount[NUM_AA_NL];
    int cNLCount[NUM_AA_NL];
    
    if (NUM_AA_NL) {
        //zero count arrays
        for (i=0; i<NUM_AA_NL; i++) {
            nNLCount[i] = 0;
            cNLCount[i] = 0;
            //initialize aa counts to y fragment
            for(iFragment=0; iFragment<cCandidate.iPeptideLength; iFragment++) {
                if (nlValuesAA[i].hasAA[cCandidate.sPeptide[iFragment]]) {
                    cNLCount[i] += 1;
                    cNLTotal += 1;
                }
            }
        }
    }
    
    // Base Masses
    nTermFragMass = GPU_MASS_PEPTIDE_NTERM;
    if (cCandidate.cBefore == '-') nTermFragMass += GPU_MASS_PROTEIN_NTERM;
    if (cCandidate.ntermMod) nTermFragMass += NTERM_MOD_MASSES[cCandidate.ntermMod-1];

    //baseMass = cCandidate.fPeptideMass - PROTON_MASS;

    //if (iThread==0)
    //   printf("%s\n", cCandidate.sPeptide);
    //loop through the peptide
    for (iFragment=0; iFragment<cCandidate.iPeptideLength-1; iFragment++) {
        aa = cCandidate.sPeptide[iFragment];
        //update phos mods
        if (NUM_AA_NL && cNLTotal) {
            for (i=0; i<NUM_AA_NL; i++) {
                if (nlValuesAA[i].hasAA[aa]) {
                    nNLCount[i] += 1;
                    cNLCount[i] -= 1;
                    nNLTotal += 1;
                    cNLTotal -= 1;
                }
            }
        }

        //calculate the fragment mass
        nTermFragMass += GPU_MASS_AA[aa];
        //baseMass -= GPU_MASS_AA[aa];

        for (iCharge=1; iCharge<iPrecursorCharge; iCharge++) {
            
            if (USE_A_IONS) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c A %f\n", aa, (nTermFragMass + PROTON_MASS*iCharge - CO_MASS) / iCharge);
            }
            if (USE_B_IONS) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c B %f\n", aa, (nTermFragMass + PROTON_MASS*iCharge) / iCharge);
            }
            if (USE_C_IONS) {
                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c C %f\n", aa, (nTermFragMass + PROTON_MASS*iCharge + NH3_MASS) / iCharge);
            }

            if (USE_X_IONS) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c X %f\n", aa, (baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS) / iCharge);
            }
            if (USE_Y_IONS) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1)) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c Y %f\n", aa, (baseMass - nTermFragMass + PROTON_MASS*(iCharge-1)) / iCharge);
            }
            if (USE_Z_IONS) {
                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                if (0 <= iBin && iBin < GPU_NUM_BINS)
                    fScore += fSpectrum[iBin];
                //if (iThread==0)
                //    printf("%c Z %f\n", aa, (baseMass -nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS) / iCharge);
            }
            
            if (NUM_AA_NL) {
                if (cNLTotal || nNLTotal){// && bMatchedInt > 0) {
                    for (i=0; i<NUM_AA_NL; i++) {
                        if (nNLCount[i]) {
                            if (USE_A_IONS) {
                                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                            if (USE_B_IONS) {
                                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                            if (USE_C_IONS) {
                                iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                        }
                        if (cNLCount[i]) {
                            if (USE_X_IONS) {
                                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                            if (USE_Y_IONS) {
                                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                            if (USE_Z_IONS) {
                                iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                                if (0 <= iBin && iBin < GPU_NUM_BINS)
                                    fScore += fSpectrum[iBin] * nlValuesAA[i].weighting;
                            }
                        }
                    }
                }
            }                     

            if (NUM_NTERM_NL){// && bMatchedInt > 0) {
                for (i=0; i<NUM_NTERM_NL; i++) {
                    if (nlValuesNterm[i].modNum != cCandidate.ntermMod)
                        continue;
                    if (USE_A_IONS) {
                        iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - CO_MASS - nlValuesNterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesNterm[i].weighting;
                    }
                    if (USE_B_IONS) {
                        iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge - nlValuesNterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesNterm[i].weighting;
                    }
                    if (USE_C_IONS) {
                        iBin = (int) ((nTermFragMass + PROTON_MASS*iCharge + NH3_MASS - nlValuesNterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesNterm[i].weighting;
                    }
                }
            }
            if (NUM_CTERM_NL){// && yMatchedInt > 0) {
                for (i=0; i<NUM_CTERM_NL; i++) {
                    if (nlValuesCterm[i].modNum != cCandidate.ctermMod)
                        continue;
                    if (USE_X_IONS) {
                        iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - H_MASS + CO_MASS - nlValuesCterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesCterm[i].weighting;
                    }
                    if (USE_Y_IONS) {
                        iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - nlValuesCterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesCterm[i].weighting;
                    }
                    if (USE_Z_IONS) {
                        iBin = (int) ((baseMass - nTermFragMass + PROTON_MASS*(iCharge-1) - NH3_MASS - nlValuesCterm[i].massDelta) / iCharge / GPU_TOLERANCE + FRAGMENT_BIN_OFFSET);
                        if (0 <= iBin && iBin < GPU_NUM_BINS)
                            fScore += fSpectrum[iBin] * nlValuesCterm[i].weighting;
                    }
                }
            }
        }
    }                               

    gpu_fScores[iThread] = fScore*PRIMARY_INTENSITY/10000.0f;
    
}

__kernel void gpu_score_reduction1(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset, __global float* fNextScore, long lNextScoreOffset, __local float* shared, __constant char* UN_MOD_AA) {
    
    int   iThread, iThreadIndex, iScoresPerThread, localSize, i;
    float f;
    __local int* bNextScore;
    __local int* iMatches;
    __local float* fScores;
    __local int* occurrences;
    __local float* nextScores;

    if (get_global_id(0) >= iNumScores)
        return;
    iThread = get_local_id(0);

    localSize = min((int)get_local_size(0), iNumScores);
    iScoresPerThread = ceil((float)iNumScores/localSize);
    iThreadIndex = iThread*iScoresPerThread;
    
    mPSMs += lPSMsOffset;
    fNextScore += lNextScoreOffset;
    
    // shared memory pointers
    iMatches =   (__local int*) shared;
    //fScores =    (float*) &iMatches[blockDim.x];
    //bNextScore = (int*)   &fScores[blockDim.x]    ;
    fScores =    (__local float*) &iMatches[localSize];
    nextScores = (__local float*) &fScores[localSize];
    occurrences = (__local int*)  &nextScores[localSize];
    bNextScore = (__local int*)   &occurrences[localSize];

    // thread, scores per thread, and thread index
    //iThread = threadIdx.x;
    //iScoresPerThread = iNumScores/blockDim.x;
    
    // Part One: Each thread reduces a portion of the scores

    // initialize local memory
    iMatches[iThread] = -1;
    fScores[iThread] = MIN_SCORE;
    nextScores[iThread] = MIN_SCORE;
    occurrences[iThread] = 1;
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // loop through scores for this thread
    for (i = iThreadIndex; i<(iThread + 1)*iScoresPerThread && (iThread + 1)*iScoresPerThread<iNumScores; i++) {
        
        //prefetch score to a register
        f = fCandidateScores[i];

        //best match
        if (f > fScores[iThread]) {
            if (GPU_FIX_DELTA_SCORE == 0 || !backbone_match(&cCandidates[iThread], &cCandidates[i]))
                nextScores[iThread] = fScores[iThread];
            iMatches[iThread] = i;
            fScores[iThread] = f;
        }
        else if (f == nextScores[iThread]) {
            if (duplicate(&cCandidates[iThread], &cCandidates[i])) {
                if (cCandidates[iThread].iProtein != cCandidates[i].iProtein) {
                    cCandidates[iThread].iProtein = min(cCandidates[iThread].iProtein, cCandidates[i].iProtein);
                    occurrences[iThread] += 1;
                }
            }
        }
        else if (f > nextScores[iThread]) {
            if (GPU_FIX_DELTA_SCORE == 0 || !backbone_match(&cCandidates[iThread], &cCandidates[i]))
                nextScores[iThread] = f;
        }    
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    for (i=localSize/2; i>=1; i/=2) {
        if (iThread<i) {
            if (fScores[iThread + i] > fScores[iThread]) {
                if (GPU_FIX_DELTA_SCORE == 0 || !backbone_match(&cCandidates[iMatches[iThread]], &cCandidates[iMatches[iThread + i]]))
                    nextScores[iThread] = fScores[iThread];
                else if (nextScores[iThread + 1] > nextScores[iThread])
                    nextScores[iThread] = nextScores[iThread + i];
                iMatches[iThread] = iMatches[iThread + i];
                fScores[iThread] = fScores[iThread + i];
            }
            else if (fScores[iThread + i] == nextScores[iThread]) {
                if (duplicate(&cCandidates[iMatches[iThread + i]], &cCandidates[iMatches[iThread]])) {
                    if (cCandidates[iMatches[iThread]].iProtein != cCandidates[iMatches[iThread + i]].iProtein) {
                        cCandidates[iMatches[iThread]].iProtein = min(cCandidates[iMatches[iThread]].iProtein, cCandidates[iMatches[iThread + i]].iProtein);
                        occurrences[iThread] += 1;
                    }
                    if (nextScores[iThread + 1] > nextScores[iThread])
                        nextScores[iThread] = nextScores[iThread + i];
                }
            }
            else if (fScores[iThread + i] > nextScores[iThread]) {
                if (GPU_FIX_DELTA_SCORE == 0 || !backbone_match(&cCandidates[iMatches[iThread]], &cCandidates[iMatches[iThread + i]]))
                    nextScores[iThread] = fScores[iThread + i];
                else if (nextScores[iThread + i] > nextScores[iThread])
                    nextScores[iThread] = nextScores[iThread + i];
            }
            else if (nextScores[iThread + 1] > nextScores[iThread])
                nextScores[iThread] = nextScores[iThread + i];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // thread 0 holds the top PSM
    if (iThread == 0) {
        *bNextScore = 0;

        //new best score
        if (fScores[0] > mPSMs[0].fScore) {
            // update next best score
            if (GPU_FIX_DELTA_SCORE == 0 || !gpu_backbone_match(&cCandidates[iMatches[0]], &mPSMs[0])) {
                *fNextScore = max(mPSMs[0].fScore, nextScores[0]);
            }

            // update top PSM
            mPSMs[0].iPeptideLength = cCandidates[iMatches[0]].iPeptideLength;
            for (i=0; i<cCandidates[iMatches[0]].iPeptideLength; i++)
                mPSMs[0].sPeptide[i] = cCandidates[iMatches[0]].sPeptide[i];
            mPSMs[0].sPeptide[i] = '\0';
            mPSMs[0].ntermMod = cCandidates[iMatches[0]].ntermMod;
            mPSMs[0].ctermMod = cCandidates[iMatches[0]].ctermMod;
            mPSMs[0].fPeptideMass = cCandidates[iMatches[0]].fPeptideMass;
            mPSMs[0].fScore = fScores[0];
            mPSMs[0].iProtein = cCandidates[iMatches[0]].iProtein;
            mPSMs[0].cBefore = cCandidates[iMatches[0]].cBefore;
            mPSMs[0].cAfter = cCandidates[iMatches[0]].cAfter;
            mPSMs[0].iNumOccurrences = occurrences[0];
        }

        //duplicate
        else if (gpu_duplicate(&cCandidates[iMatches[0]], fScores[0], &mPSMs[0])) {
            if (cCandidates[iMatches[0]].iProtein != mPSMs[0].iProtein) {
                mPSMs[0].iProtein = min(mPSMs[0].iProtein, cCandidates[iMatches[0]].iProtein);
                mPSMs[0].iNumOccurrences += occurrences[0];
            }
            if (nextScores[0] > *fNextScore)
                *fNextScore = nextScores[0];
        }

        //next best score
        else if (fScores[0] > *fNextScore)
                *fNextScore = fScores[0];
        else if (nextScores[0] > *fNextScore)
                *fNextScore = nextScores[0];

        fCandidateScores[iMatches[0]] = MIN_SCORE;
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);

    /*
    while (*bNextScore == 0) {
    
        // Part One: Each thread reduces a portion of the scores

        // initialize best match
        iMatches[iThread] = -1;
        fScores[iThread] = 0.0f;

        for (i = iThreadIndex; i < (iThread + 1) * iScoresPerThread; i++) {
            //prefetch
            f = fCandidateScores[i];

            //update best match
            if (f > fScores[iThread]) {
                iMatches[iThread] = i;
                fScores[iThread] = f;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        for (i=BLOCKDIM_REDUCE/2; i>=1; i/=2) {
            if (iThread<i)
                gpu_reduce_score(fScores, iMatches, i);
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        
        // thread 0 holds the next best PSM
        if (iThread == 0) {
            // keep the current next best score
            if (fScores[0] <= *fNextScore) {
                *bNextScore = 1;
            }

            // duplicate
            else if (gpu_duplicate(&cCandidates[iMatches[0]], fScores[0], &mPSMs[0])) {
                if (cCandidates[iMatches[0]].iProtein != mPSMs[0].iProtein) mPSMs[0].iNumOccurrences += 1;
            }

            // new next best score
            else if (GPU_FIX_DELTA_SCORE == 0 || gpu_backbone_match(&cCandidates[iMatches[0]], &mPSMs[0], UN_MOD_AA) == 0) {
                *fNextScore = fScores[0];
                *bNextScore = 1;
            }

            fCandidateScores[iMatches[0]] = MIN_SCORE;
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
    */
}
/*
 *  Score reduction kernel (based on cuda sdk kernel6)
 */

__kernel void gpu_score_reduction2(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset, __global float* fNextScore, long lNextScoreOffset, __local float* shared, __constant char* UN_MOD_AA)
{
    int   iThread, iThreadIndex, iScoresPerThread, i, j, k, c;
    float score;
    
    iThread = get_local_id(0);
    if (iThread > 0) return;

    __local int topNMatches[10];
    __local float topNScores[10];
    __local float minScore;
    minScore = MIN_SCORE;

    for (i=0; i<10; i++) {
        topNScores[i] = MIN_SCORE;
    }
    
    for (i=0; i<4096; i++) {
        //prefetch score to a register
        score = fCandidateScores[i];
        if (score < minScore)
            continue;
        
        for (j=0; j<10; j++) {
            if (score > topNScores[j]) {
                //shift right and insert
                for (k=10-2; k>=j; k--) {
                    topNMatches[k+1] = topNMatches[k];
                    topNScores[k+1] = topNScores[k];
                }
                topNMatches[j] = i;
                topNScores[j] = score;
                minScore = topNScores[10-1];
                break;
            }
        }
    }

    for (i=0; i<10; i++) {
        for (j=0; j<10; j++) {
            if (topNScores[i] > mPSMs[j].fScore) {
                //shift right and insert
                for (k=10-2; k>=j; k--) {
                    mPSMs[k+1] = mPSMs[k];
                }
                mPSMs[j].iPeptideLength = cCandidates[topNMatches[i]].iPeptideLength;
                for (int c=0; c<cCandidates[topNMatches[i]].iPeptideLength; c++)
                    mPSMs[i].sPeptide[c] = cCandidates[topNMatches[i]].sPeptide[c];
                mPSMs[j].sPeptide[i] = '\0';
                mPSMs[j].ntermMod = cCandidates[topNMatches[i]].ntermMod;
                mPSMs[j].ctermMod = cCandidates[topNMatches[i]].ctermMod;
                mPSMs[j].fPeptideMass = cCandidates[topNMatches[i]].fPeptideMass;
                mPSMs[j].fScore = topNScores[i];
                mPSMs[j].iProtein = cCandidates[topNMatches[i]].iProtein;
                mPSMs[j].cBefore = cCandidates[topNMatches[i]].cBefore;
                mPSMs[j].cAfter = cCandidates[topNMatches[i]].cAfter;
                mPSMs[j].iNumOccurrences = 1;
                break;
            }
        }
    }

    fNextScore[lNextScoreOffset] = mPSMs[1].fScore;
}

__kernel void gpu_score_reduction(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset, long lNextScoreOffset, __local float* shared)
{
    int   i, j;
    float f;
    int iThread = get_local_id(0);
    int localSize = get_local_size(0);
    __local int* iMatches;
    //__local int* occurrences;
    __local float* fScores;
    //__local mObj psms[NUM_OUTPUT_PSMS];
    __local int psmInd;
    //__local float currentScore;
    psmInd = 0;
    
    mPSMs += lPSMsOffset;

    //currentScore = mPSMs[NUM_OUTPUT_PSMS-1].fScore;
    //copyEvent = async_work_group_copy((__local char*)psms, (__global char*)mPSMs, NUM_OUTPUT_PSMS*80, 0);
    //wait_group_events(1, &copyEvent);
    
    // shared memory pointers
    iMatches =   (__local int*) shared;
    fScores =    (__local float*) (iMatches + get_local_size(0));

    //__local int topNMatches[NUM_OUTPUT_PSMS];
    //__local float topNScores[NUM_OUTPUT_PSMS];
    //float lowestScore = mPSMs[9].fScore;

    // thread, scores per thread, and thread index
    int iScoresPerThread = iNumScores/localSize;
    int iThreadIndex = iThread*iScoresPerThread;
    
    // Part One: Each thread reduces a portion of the scores
    while (psmInd<NUM_OUTPUT_PSMS) {
        // initialize local memory
        iMatches[iThread] = iThreadIndex;
        fScores[iThread] = fCandidateScores[iThreadIndex];

        // loop through scores for this thread
        for (i = iThreadIndex+1; i < (iThread + 1) * iScoresPerThread; i++) {
            //prefetch score to a register
            f = fCandidateScores[i];
            //best match
            if (f > fScores[iThread]) {
                iMatches[iThread] = i;
                fScores[iThread] = f;
            }
            // else if (GPU_TRACK_DUPLICATES && f == fScores[iThread]) {
            //     if (duplicate(
                
                
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        for (i=localSize/2; i>=1; i/=2) {
            if (iThread<i) {
                //gpu_reduce_score(fScores, iMatches, iThread, i);
                if (fScores[iThread] < fScores[iThread + i]) {
                    iMatches[iThread] = iMatches[iThread + i];
                    fScores[iThread] = fScores[iThread + i];
                }
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // thread 0 holds the top PSM
        if (iThread == 0) {
            if (fScores[0] < mPSMs[NUM_OUTPUT_PSMS-1].fScore)
                psmInd = NUM_OUTPUT_PSMS; //break from loops
            
            //new best score
            for (; psmInd<NUM_OUTPUT_PSMS; psmInd++) {
                if (fScores[0] > mPSMs[psmInd].fScore) {
                    // update next best score
                    //if (GPU_FIX_DELTA_SCORE == 0 || gpu_backbone_match(&cCandidates[iMatches[0]], &mPSMs[0], UN_MOD_AA) == 0) {
                    //    *fNextScore = mPSMs[0].fScore;
                    //}

                    //shift right
                    for (j=NUM_OUTPUT_PSMS-1; j>psmInd; j--)
                        mPSMs[j] = mPSMs[j-1];

                    // update ith PSM
                    mPSMs[psmInd].iPeptideLength = cCandidates[iMatches[0]].iPeptideLength;
                    int c = 0;
                    for (c=0; c<cCandidates[iMatches[0]].iPeptideLength; c++)
                        mPSMs[psmInd].sPeptide[c] = cCandidates[iMatches[0]].sPeptide[c];
                    mPSMs[psmInd].sPeptide[c] = '\0';
                    mPSMs[psmInd].ntermMod = cCandidates[iMatches[0]].ntermMod;
                    mPSMs[psmInd].ctermMod = cCandidates[iMatches[0]].ctermMod;
                    mPSMs[psmInd].fPeptideMass = cCandidates[iMatches[0]].fPeptideMass;
                    mPSMs[psmInd].fScore = fScores[0];
                    mPSMs[psmInd].iProtein = cCandidates[iMatches[0]].iProtein;
                    mPSMs[psmInd].cBefore = cCandidates[iMatches[0]].cBefore;
                    mPSMs[psmInd].cAfter = cCandidates[iMatches[0]].cAfter;
                    mPSMs[psmInd].iNumOccurrences = 1;

                    //currentScore = fScores[0];
                    break;
                }
                else if (gpu_duplicate(&cCandidates[iMatches[0]], fScores[0], &mPSMs[psmInd]))
                    break;

                //duplicate
                //else if (gpu_duplicate(&cCandidates[iMatches[0]], fScores[0], &mPSMs[0])) {
                //    if (cCandidates[iMatches[0]].iProtein != mPSMs[0].iProtein) mPSMs[0].iNumOccurrences += 1;
                //}

                //next best score
                //else {
                //    if (fScores[0] > *fNextScore) {
                //        *fNextScore = fScores[0];
                //    }
            
                //    *bNextScore = 1;
                //}
            }
            fCandidateScores[iMatches[0]] = MIN_SCORE;
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //copyEvent = async_work_group_copy((__global char*)mPSMs, (__local char*)psms, NUM_OUTPUT_PSMS*80, 0);
    //wait_group_events(1, &copyEvent);
    //if (iThread == 0)
    //    printf("iteration=%d\n", iteration);
}


//=================================================================================================
// Reduction Kernel Functions (called from global kernels)
//=================================================================================================

void gpu_reduce_score(__local float* fScores, __local int* iMatches, int iThread, int iThreadOffset)
{
    if (fScores[iThread] < fScores[iThread + iThreadOffset]) {
        iMatches[iThread] = iMatches[iThread + iThreadOffset];
        fScores[iThread] = fScores[iThread + iThreadOffset];
    }
}

bool gpu_duplicate(__global cObj *cCandidate, float fCandidateScore, __global mObj* mPSM)
{
    // compare the score, length, and string
    if (fCandidateScore != mPSM->fScore) return 0;
    if (cCandidate->iPeptideLength != mPSM->iPeptideLength) return 0;
    return gpu_strcmp(cCandidate->sPeptide, mPSM->sPeptide, cCandidate->iPeptideLength);
}

bool duplicate(__global cObj* cCandidate1, __global cObj* cCandidate2)
{
    // compare the length and string
    if (cCandidate1->iPeptideLength != cCandidate2->iPeptideLength) return 0;
    return gpu_strcmp(cCandidate1->sPeptide, cCandidate2->sPeptide, cCandidate1->iPeptideLength);
}

bool gpu_backbone_match(__global cObj *cCandidate, __global mObj* mPSM)
{
    if (cCandidate->iPeptideLength != mPSM->iPeptideLength) return 0;
    return gpu_strcmpi(cCandidate->sPeptide, mPSM->sPeptide, cCandidate->iPeptideLength);
}

bool backbone_match(__global cObj* cCandidate1, __global cObj* cCandidate2) {
    if (cCandidate1->iPeptideLength != cCandidate2->iPeptideLength) return 0;
    return gpu_strcmpi(cCandidate1->sPeptide, cCandidate2->sPeptide, cCandidate1->iPeptideLength);
}

//=================================================================================================
// Kernel Utility Functions
//=================================================================================================

//normal strcmp
bool gpu_strcmp(__global const unsigned char* s1, __global const unsigned char* s2, int peptideLength)
{
    int i=0;
    while (i<peptideLength && (s1[i] == s2[i]))
        i++;
    return s1[i] == s2[i];
}

//case insensitive strcmp - only works for purely alphabetic strings
//int gpu_strcmpi(__global const unsigned char * s1, __global const unsigned char * s2, __constant char* UN_MOD_AA)
// {
//     for (; *s1 == *s2 || *s1+CASE_SHIFT == *s2 || *s1 == *s2+CASE_SHIFT; ++s1, ++s2)
//         if (*s1 == 0) return 0;
//     return *s1 - *s2;
//}
bool gpu_strcmpi2(__global const unsigned char* s1, __global const unsigned char* s2)
 {
     while (*s1 && (((*s1)&31) == ((*s2)&31))) {
         ++s1;
         ++s2;
     }
     return ((*s1)&31) == ((*s2)&31);
 }

bool gpu_strcmpi(__global const unsigned char* s1, __global const unsigned char* s2, int peptideLength)
 {
     int i=0;
     while (i<peptideLength && (((s1[i]%32) == (s2[i]%32))))
         i++;
     return (s1[i]%32) == (s2[i]%32);
 }

__kernel void cl_memset(__global char * dst, int c, unsigned long n)
{
    int global_id = get_global_id(0);
    if(global_id < n)
        dst[global_id] = c;
}
