/*
 *  kernels.cl
 *
 *  Created by Brendan Faherty on 8/31/10.
 *  Copyright 2010 Dartmouth College. All rights reserved.
 *
 */

#define STRING_SIZE 256
#define CASE_SHIFT 32
#define aa2i(c) (int) (c - 'A')

//MEA: include tempest.h instead of copying these?
// scoring
#define PRIMARY_INTENSITY 50
#define PHOS_NL_INTENSITY 25
#define FLANKING_INTENSITY 25
#define NL_INTENSITY 10
#define MIN_SCORE -1000.0f

// constants
#define PROTON_MASS     1.00727646688
#define H_MASS          1.007825032
#define C_MASS          12.0
#define N_MASS          14.003074007
#define O_MASS          15.994914622

#define OH_MASS         17.002739654
#define NH3_MASS        17.026549103
#define H2O_MASS        18.010564686
#define CO_MASS         27.994914622
#define H3PO4_MASS      97.97689509

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

//Mass delta and score weighting for a neutral loss event
typedef struct nlValue {
    int   numSites;
    char  aaString[64];
    float massDelta;
    float weighting;
} nlValue;

void gpu_reduce_score(__local float* fScores, __local int* iMatches, int iThreadOffset);
bool gpu_duplicate(__global cObj *cCandidate, float fCandidateScore, __global mObj* mPSM);
bool gpu_backbone_match(__global cObj *cCandidate, __global mObj* mPSM);
int gpu_strcmp(__global const char * s1, __global const char * s2);
int gpu_strcmpi(__global const char * s1, __global const char * s2);

//=================================================================================================
// Global Kernel Functions
//=================================================================================================


__kernel void gpu_build(__global float* fSpectrum, long lSpectrumOffset, int iPeakCount, __global int* iPeakBins, __global float* fPeakInts, long lPeakOffset)
{
    //printf("gpu_build\n");
    int iBin, iThread;
    float intensity;
    //iThread = blockIdx.x * blockDim.x + threadIdx.x;
    iThread = get_global_id(0);
    
    fSpectrum += lSpectrumOffset;
    iPeakBins += lPeakOffset;
    fPeakInts += lPeakOffset;

    if (iThread < iPeakCount) {
        iBin = iPeakBins[iThread];
        intensity = fPeakInts[iThread];
        fSpectrum[iBin] = intensity;
    }

#if GPU_FLANKING
    barrier(CLK_GLOBAL_MEM_FENCE);
    if (iThread < iPeakCount) {
        if (iBin > 0)
            intensity += fSpectrum[iBin-1]/2;
        if (iBin < GPU_NUM_BINS-1)
            intensity += fSpectrum[iBin+1]/2;
    }
    barrier(CLK_GLOBAL_MEM_FENCE);

    if (iThread < iPeakCount)
        fSpectrum[iBin] = intensity;
#endif
    
}

/*
 * Kernel for transforming a single observed spectrum in global memory.
 * 
 *      fScratch        Built (rasterized) spectrum, prior to transformation.
 *      fSpectrum       Memory for transformed spectrum.
 */

__kernel void gpu_transform_scratch(__global float* fScratch, long lScratchOffset, __global float* fSpectrum, long lSpectrumOffset)
{

    int iStart = get_global_id(0);
    if (iStart >= GPU_NUM_BINS) return;

    fScratch += lScratchOffset;
    fSpectrum += lSpectrumOffset;

    float fIntensity = 0.0f;

    // Calculate the initial fIntensity for each thread
    for (int iBin = iStart-GPU_XWIDTH; iBin <= iStart+GPU_XWIDTH; iBin++) {
        fIntensity += fScratch[iBin];
    }
    
    fSpectrum[iStart] = (float) fScratch[iStart] - (fIntensity)/(GPU_XWIDTH*2+1);
        
    //update the intensity
    //fIntensity -= fScratch[iStart-GPU_XWIDTH];
    //fIntensity += fScratch[iStart+GPU_XWIDTH+1];
}

__kernel void gpu_transform(__global float* fSpectrum, long lSpectrumOffset)
{

    int iStart = get_global_id(0);
    if (iStart >= GPU_NUM_BINS) return;

    fSpectrum += lSpectrumOffset;

    float fIntensity = 0.0f;

    // Calculate the initial fIntensity for each thread
    int windowStart = max(iStart-GPU_XWIDTH, 0);
    int windowEnd = min(iStart+GPU_XWIDTH, GPU_NUM_BINS-1);
    for (int iBin = windowStart; iBin <= windowEnd; iBin++) {
        fIntensity += fSpectrum[iBin];
    }
    
    fSpectrum[iStart] -= (fIntensity)/(GPU_XWIDTH*2+1);
        
    //update the intensity
    //fIntensity -= fScratch[iStart-GPU_XWIDTH];
    //fIntensity += fScratch[iStart+GPU_XWIDTH+1];
}

__kernel void gpu_transform_local(__global float* fSpectrum, long lSpectrumOffset)
{

    int iStart = get_group_id(0);
    if (iStart >= GPU_NUM_BINS) return;
    int localIndex = iStart + (get_local_id(0) - GPU_XWIDTH/2);
    if (localIndex >= iStart)
        localIndex += 1;

    fSpectrum += lSpectrumOffset;

    __local float fIntensity;
    fIntensity = 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (localIndex >= 0 && localIndex < GPU_NUM_BINS)
        atomic_xchg(&fIntensity, fIntensity+fSpectrum[localIndex]);

    barrier(CLK_LOCAL_MEM_FENCE);
    barrier(CLK_GLOBAL_MEM_FENCE);
          
    fSpectrum[iStart] -= (fIntensity)/(GPU_XWIDTH*2+1);
        
    //update the intensity
    //fIntensity -= fScratch[iStart-GPU_XWIDTH];
    //fIntensity += fScratch[iStart+GPU_XWIDTH+1];
}

__kernel void gpu_transform2(__global float* fScratch, long lScratchOffset, __global float* fSpectrum, long lSpectrumOffset)
{
    int iBin, iBinsPerThread, iStart;
    float fIntensity;

    //if (blockIdx.x > 0) return;
    if (get_group_id(0) > 0) return;

    fScratch += lScratchOffset;
    fSpectrum += lSpectrumOffset;

    fIntensity = 0.0f;

    //iBinsPerThread = GPU_NUM_BINS / blockDim.x;
    //iStart = threadIdx.x * iBinsPerThread;
    iBinsPerThread = GPU_NUM_BINS / get_local_size(0);
    iStart = get_local_id(0) * iBinsPerThread;

    // Calculate the initial fIntensity for each thread
    for (iBin = iStart-GPU_XWIDTH; iBin <= iStart+GPU_XWIDTH; iBin++) {
        fIntensity += fScratch[iBin];
    }
    
    // Loop through the bins for each thread
    for (iBin=iStart; iBin+1<iStart+iBinsPerThread; iBin++) {
        //store the new intensity
        fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(GPU_XWIDTH*2+1);
        
        //update the intensity
        fIntensity -= fScratch[iBin-GPU_XWIDTH];
        fIntensity += fScratch[iBin+GPU_XWIDTH+1];
    }
    
    // use the last thread (if necessary) for the rest of the bins
    //if (threadIdx.x == blockDim.x-1 && fIntensity>0.0f) {
    if (get_local_id(0) == get_local_size(0)-1 && fIntensity>0.0f) {
        for (; iBin+1<GPU_NUM_BINS; iBin++) {
            //store the new intensity
            fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(GPU_XWIDTH*2+1);
            
            //update the intensity
            fIntensity -= fScratch[iBin-GPU_XWIDTH];
            fIntensity += fScratch[iBin+GPU_XWIDTH+1];
        }
    }
    
    // and get the last one
    fSpectrum[iBin] = (float) fScratch[iBin] - (fIntensity)/(GPU_XWIDTH*2+1);
}

__kernel void gpu_score(int iPrecursorCharge, int iNumCandidates, __global cObj* gpu_cCandidates, __global float* gpu_fScores, __global float* gpu_fSpectra, long lSpectrumOffset,
                        __constant float* GPU_MASS_AA, __constant float* GPU_MASS_PROTON, __constant nlValue* nlValuesNterm, __constant nlValue* nlValuesCterm, __constant nlValue* nlValuesAA)
{
    //printf("gpu_score\n");
    int iThread = get_global_id(0);
    //iThread = blockIdx.x*blockDim.x + threadIdx.x;
    if (iThread >= iNumCandidates) {
        gpu_fScores[iThread] = 0.0f;
        return;
    }
  
    __global float* fSpectrum = gpu_fSpectra + lSpectrumOffset;
    __local cObj cCandidates[BLOCKDIM_SCORE];
    cCandidates[get_local_id(0)] = gpu_cCandidates[get_group_id(0)*BLOCKDIM_SCORE+get_local_id(0)];

    __local cObj* cCandidate = &(cCandidates[get_local_id(0)]);

    int iCharge;
    int iBin;
    int iFragment;
    float fBMass;
    float fYMassBase;
    float fScoreP = 0.0f;
    float fScoreNL = 0.0f;
    int bNLTotal = 0;
    int yNLTotal = 0;
    int bNLCount[NUM_AA_NL];
    int yNLCount[NUM_AA_NL];
    float bMatchedInt;
    float yMatchedInt;

    if (NUM_AA_NL) {
        //zero count arrays
        for (int i=0; i<NUM_AA_NL; i++) {
            bNLCount[i] = 0;
            yNLCount[i] = 0;
            //initialize aa counts to y fragment
            for(iFragment=0; iFragment<cCandidate->iPeptideLength; iFragment++) {
                char aa = cCandidate->sPeptide[iFragment];
                for (int j=0; j<nlValuesAA[i].numSites; j++) {
                    if (nlValuesAA[i].aaString[j] == aa) {
                        yNLCount[i] += 1;
                        yNLTotal += 1;
                        break;
                    }
                }
            }
        }
    }
    
    // Base Masses
    fBMass = GPU_MASS_PEPTIDE_NTERM;
    if (cCandidate->cBefore == '-') fBMass += GPU_MASS_PROTEIN_NTERM;
    if (cCandidate->bNtermMod) fBMass += GPU_MASS_NTERM_MOD;

    fYMassBase = cCandidate->fPeptideMass - PROTON_MASS;

    //loop through the peptide
    for (iFragment=0; iFragment<cCandidate->iPeptideLength-1; iFragment++) {
        //update phos mods
        if (NUM_AA_NL && yNLTotal) {
            char aa = cCandidate->sPeptide[iFragment];
            for (int i=0; i<NUM_AA_NL; i++) {
                for (int j=0; j<nlValuesAA[i].numSites; j++) {
                    if (nlValuesAA[i].aaString[j] == aa) {
                        bNLCount[i] += 1;
                        yNLCount[i] -= 1;
                        bNLTotal += 1;
                        yNLTotal -= 1;
                        break;
                    }
                }
            }
        }

        //calculate the fragment mass
        fBMass += GPU_MASS_AA[aa2i(cCandidate->sPeptide[iFragment])];

        for (iCharge=1; iCharge<iPrecursorCharge; iCharge++) {
            // B Primary
            iBin = (int) ((fBMass + GPU_MASS_PROTON[iCharge]) / iCharge / GPU_TOLERANCE + 0.5);
            if (0 <= iBin && iBin < GPU_NUM_BINS) {
                bMatchedInt = fSpectrum[iBin];
                fScoreP += bMatchedInt;
            }
            
            // Y Primary
            iBin = (int) ((fYMassBase - fBMass + GPU_MASS_PROTON[iCharge]) / iCharge / GPU_TOLERANCE + 0.5);
            if (0 <= iBin && iBin < GPU_NUM_BINS) {
                yMatchedInt = fSpectrum[iBin];
                fScoreP += yMatchedInt;
            }

            if (NUM_AA_NL) {
                // NLs from B ion
	      if (bNLTotal){// && bMatchedInt > 0) {
                    for (int i=0; i<NUM_AA_NL; i++) {
                        if (bNLCount[i]) {
                            iBin = (int) ((fBMass + GPU_MASS_PROTON[iCharge] - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + 0.5f);
                            if (0 <= iBin && iBin < GPU_NUM_BINS)
                                fScoreNL += fSpectrum[iBin] * nlValuesAA[i].weighting;
                        }
                    }
                }
	      if (yNLTotal){// && yMatchedInt > 0) {
                    for (int i=0; i<NUM_AA_NL; i++) {
                        if (yNLCount[i]) {
                            iBin = (int) ((fYMassBase - fBMass + GPU_MASS_PROTON[iCharge] - nlValuesAA[i].massDelta) / iCharge / GPU_TOLERANCE + 0.5);
                            if (0 <= iBin && iBin < GPU_NUM_BINS) 
                                fScoreNL += fSpectrum[iBin] * nlValuesAA[i].weighting;
                        }
                    }
                }
            }                     

            if (NUM_NTERM_NL){// && bMatchedInt > 0) {
                for (int i=0; i<NUM_NTERM_NL; i++) {
                    iBin = (int) ((fBMass + GPU_MASS_PROTON[iCharge] - nlValuesNterm[i].massDelta) / iCharge / GPU_TOLERANCE + 0.5f);
                    if (0 <= iBin && iBin < GPU_NUM_BINS)
                        fScoreNL += fSpectrum[iBin] * nlValuesNterm[i].weighting;
                }
            }
            if (NUM_CTERM_NL){// && yMatchedInt > 0) {
                for (int i=0; i<NUM_CTERM_NL; i++) {
                    iBin = (int) ((fYMassBase - fBMass + GPU_MASS_PROTON[iCharge] - nlValuesCterm[i].massDelta) / iCharge / GPU_TOLERANCE + 0.5);
                    if (0 <= iBin && iBin < GPU_NUM_BINS)
                        fScoreNL += fSpectrum[iBin] * nlValuesCterm[i].weighting;
                }
            }
        }
    }                                 
    
    gpu_fScores[iThread] = (fScoreP + fScoreNL)*PRIMARY_INTENSITY/10000.0f;
    
}

/*
 *  Score reduction kernel (based on cuda sdk kernel6)
 */

__kernel void gpu_score_reduction(int iNumScores, __global cObj* cCandidates, __global float* fCandidateScores, __global mObj* mPSMs, long lPSMsOffset, __global float* fNextScore, long lNextScoreOffset, __local float* shared)
{
    int   iThread, iThreadIndex, iScoresPerThread, i;
    float f;
    __local int *bNextScore;
    __local int *iMatches;
    __local float *fScores;
    
    mPSMs += lPSMsOffset;
    fNextScore += lNextScoreOffset;
    
    // shared memory pointers
    iMatches =   (__local int*) shared;
    //fScores =    (float*) &iMatches[blockDim.x];
    //bNextScore = (int*)   &fScores[blockDim.x]    ;
    fScores =    (__local float*) &iMatches[get_local_size(0)];
    bNextScore = (__local int*)   &fScores[get_local_size(0)];

    // thread, scores per thread, and thread index
    //iThread = threadIdx.x;
    //iScoresPerThread = iNumScores/blockDim.x;
    iThread = get_local_id(0);
    iScoresPerThread = iNumScores/get_local_size(0);
    iThreadIndex = iThread*iScoresPerThread;
    
    // Part One: Each thread reduces a portion of the scores

    // initialize local memory
    iMatches[iThread] = -1;
    fScores[iThread] = MIN_SCORE;
    if (iThread == 0)
        *bNextScore = 0;
    barrier(CLK_LOCAL_MEM_FENCE);
    

    // loop through scores for this thread
    for (i = iThreadIndex; i < (iThread + 1) * iScoresPerThread; i++) {
        //prefetch score to a register
        f = fCandidateScores[i];
        
        //best match
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

    // thread 0 holds the top PSM
    if (iThread == 0) {
        *bNextScore = 0;

        //new best score
        if (fScores[0] > mPSMs[0].fScore) {
            // update next best score
            if (GPU_FIX_DELTA_SCORE == 0 || gpu_backbone_match(&cCandidates[iMatches[0]], &mPSMs[0]) == 0) {
                *fNextScore = mPSMs[0].fScore;
            }

            // update top PSM
            mPSMs[0].iPeptideLength = cCandidates[iMatches[0]].iPeptideLength;
            for (i=0; i<cCandidates[iMatches[0]].iPeptideLength; i++)
                mPSMs[0].sPeptide[i] = cCandidates[iMatches[0]].sPeptide[i];
            mPSMs[0].sPeptide[i] = '\0';
            mPSMs[0].bNtermMod = cCandidates[iMatches[0]].bNtermMod;
            mPSMs[0].bCtermMod = cCandidates[iMatches[0]].bCtermMod;
            mPSMs[0].fPeptideMass = cCandidates[iMatches[0]].fPeptideMass;
            mPSMs[0].fScore = fScores[0];
            mPSMs[0].iProtein = cCandidates[iMatches[0]].iProtein;
            mPSMs[0].cBefore = cCandidates[iMatches[0]].cBefore;
            mPSMs[0].cAfter = cCandidates[iMatches[0]].cAfter;
            mPSMs[0].iNumOccurrences = 1;
        }

        //duplicate
        else if (gpu_duplicate(&cCandidates[iMatches[0]], fScores[0], &mPSMs[0])) {
            if (cCandidates[iMatches[0]].iProtein != mPSMs[0].iProtein) mPSMs[0].iNumOccurrences += 1;
        }

        //next best score
        else {
            if (fScores[0] > *fNextScore) {
                *fNextScore = fScores[0];
            }
            
            *bNextScore = 1;
        }

        fCandidateScores[iMatches[0]] = MIN_SCORE;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

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
            else if (GPU_FIX_DELTA_SCORE == 0 || gpu_backbone_match(&cCandidates[iMatches[0]], &mPSMs[0]) == 0) {
                *fNextScore = fScores[0];
                *bNextScore = 1;
            }

            fCandidateScores[iMatches[0]] = MIN_SCORE;
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
}

//=================================================================================================
// Reduction Kernel Functions (called from global kernels)
//=================================================================================================

void gpu_reduce_score(__local float* fScores, __local int* iMatches, int iThreadOffset)
{
    if (fScores[get_local_id(0)] < fScores[get_local_id(0) + iThreadOffset]) {
        iMatches[get_local_id(0)] = iMatches[get_local_id(0) + iThreadOffset];
        fScores[get_local_id(0)] = fScores[get_local_id(0) + iThreadOffset];
    }
}

bool gpu_duplicate(__global cObj *cCandidate, float fCandidateScore, __global mObj* mPSM)
{
    // compare the score, length, and string
    if (fCandidateScore != mPSM->fScore) return 0;
    if (cCandidate->iPeptideLength != mPSM->iPeptideLength) return 0;
    if (gpu_strcmp(cCandidate->sPeptide, mPSM->sPeptide) != 0) return 0;
    return 1;
}

bool gpu_backbone_match(__global cObj *cCandidate, __global mObj* mPSM)
{
    if (cCandidate->iPeptideLength != mPSM->iPeptideLength) return 0;
    if (gpu_strcmpi(cCandidate->sPeptide, mPSM->sPeptide) != 0) return 0;
    return 1;
}


//=================================================================================================
// Kernel Utility Functions
//=================================================================================================

//normal strcmp
int gpu_strcmp(__global const char * s1, __global const char * s2)
{
    for (; *s1 == *s2; ++s1, ++s2)
        if (*s1 == 0) return 0;
    return *s1 - *s2;
}

//case insensitive strcmp - only works for purely alphabetic strings
int gpu_strcmpi(__global const char * s1, __global const char * s2)
{
    for (; *s1 == *s2 || *s1+CASE_SHIFT == *s2 || *s1 == *s2+CASE_SHIFT; ++s1, ++s2)
        if (*s1 == 0) return 0;
    return *s1 - *s2;
}

__kernel void cl_memset(__global char * dst, int c, unsigned long n)
{
    int global_id = get_global_id(0);
    if(global_id < n)
        dst[global_id] = c;
}

