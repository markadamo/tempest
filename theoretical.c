/*
 *  theoretical.c
 *
 *  Created by Brendan Faherty on 12/8/08.
 *  Copyright 2008 Dartmouth College. All rights reserved.
 *
 */

#include "tempest.h"

#define LENGTH_INPUT_BUFFER 1024

//functions
int count_references(FILE*);
void parse_reference(char*, char*);
void digest_protein(int, int, char*);
void gen_candidates(cObj, bool, bool, int*, int, int);
void store_candidate(cObj cCandidate);

/*
 *  Parse fasta database.
 */

extern void search_fasta_database()
{
    int i;
    int iProtein;
    int iLengthProtein;
    int iLengthBuffer;
    int iProteinBuffer = INITIAL_PROTEIN_BUFFER;
    char sLineBuffer[LENGTH_INPUT_BUFFER+1];
    char sProteinReference[MAX_LENGTH_REFERENCE+1];
    char *sProteinSequence;
    FILE *fp;
    eObj *e;
    float fMaxModMass;
    float fMaxMassDiff;
    float fProgress;

    if (args.bNoDigest) return;

    // Setup

    // allocate working memory
    if (0 == (sProteinSequence = (char*) malloc(iProteinBuffer*sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for the active protein sequence.\n");
        tempest_exit(EXIT_FAILURE);
    }

    // calculate peptide mass range
    fMaxModMass = 0.0f;
    for (i=0; i<params.iNumMods; i++) {
        if (params.tMods[i].cSymbol && params.tMods[i].dMassDiff > fMaxModMass)
            fMaxModMass = params.tMods[i].dMassDiff;
    }

    if (params.iModificationsMax == 0) {
        fMaxMassDiff = 0.0f;
    }
    else if (params.iModificationsMax == 1) {
        fMaxMassDiff = fMaxModMass;
        if (params.dVariableNtermMassDiff > fMaxMassDiff) fMaxMassDiff = params.dVariableNtermMassDiff;
        if (params.dVariableCtermMassDiff > fMaxMassDiff) fMaxMassDiff = params.dVariableCtermMassDiff;
    }
    else {
        fMaxMassDiff = params.iModificationsMax*fMaxModMass;
        if (params.dVariableNtermMassDiff > fMaxMassDiff) fMaxMassDiff += params.dVariableNtermMassDiff - fMaxModMass;
        if (params.dVariableCtermMassDiff > fMaxMassDiff) fMaxMassDiff += params.dVariableCtermMassDiff - fMaxModMass;
    }

    if (params.bPrecursorTolerancePPM) {
        tempest.dMinPeptideMass = tempest.fMinPrecursorMass - (tempest.fMinPrecursorMass*params.fPrecursorTolerance/1000000.0f) - fMaxMassDiff;
        tempest.dMaxPeptideMass = tempest.fMaxPrecursorMass + (tempest.fMaxPrecursorMass*params.fPrecursorTolerance/1000000.0f);
    }
    else {
        tempest.dMinPeptideMass = tempest.fMinPrecursorMass - params.fPrecursorTolerance - fMaxMassDiff;
        tempest.dMaxPeptideMass = tempest.fMaxPrecursorMass + params.fPrecursorTolerance;
    }
    
    // open fasta file
    if (0 == (fp = (FILE *) fopen(args.sFasta, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open fasta database: %s\n", args.sFasta);
        tempest_exit(EXIT_FAILURE);
    }

    // count proteins
    tempest.iNumProteins = count_references(fp);
    
    if (tempest.iNumProteins == 0) {
        fprintf(stderr, "\nERROR\tNo proteins found in fasta database: %s\n", args.sFasta);
        tempest_exit(EXIT_FAILURE);
    }

    // progress
    if (args.iPrintLevel && !args.bPrintCandidates) {
        printf(" » Digesting %d proteins...     ", tempest.iNumProteins);
    }
    
    // Get memory for the protein references
    if (0 == (sProteinReferences = (char *) malloc(tempest.iNumProteins * MAX_LENGTH_REFERENCE * sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for protein references.\n");
        tempest_exit(EXIT_FAILURE);
    }
    
    // parse fasta line by line
    iProtein = 0;
    iLengthProtein = 0;
    sProteinSequence[0] = '\0';
    while (fgets(sLineBuffer, LENGTH_INPUT_BUFFER, fp)) {
        if (strlen(sLineBuffer) <= 0) continue;
        
        // reference line
        if ('>' == sLineBuffer[0]) {
            // check for previous protein
            if (iLengthProtein > 0) {
                // process
                strcpy(&sProteinReferences[iProtein * MAX_LENGTH_REFERENCE * sizeof(char)], sProteinReference);
                digest_protein(iProtein++, iLengthProtein, sProteinSequence);

                if (args.iPrintLevel && !args.bPrintCandidates) {
                    fProgress = (100.0f * iProtein) / tempest.iNumProteins;
                    if (fProgress - floor(fProgress) < 0.01f) printf("\b\b\b\b%*d%%", 3, (int) fProgress);
                    fflush(0);
                }
            
                // reset
                iLengthProtein = 0;
                sProteinSequence[0] = '\0';
            }

            // get next reference
            parse_reference(sLineBuffer, sProteinReference);
        }

        // protein sequence line
        else {
            //get length
            iLengthBuffer = (int) strlen(sLineBuffer);
            
            //strip newlines
            while(iLengthBuffer>0 && !isalnum(sLineBuffer[iLengthBuffer-1])) {
                iLengthBuffer--;
            }

            if (iLengthBuffer <= 0)
                continue;

            // check new sequence length against buffer size
            if (iLengthProtein + iLengthBuffer > iProteinBuffer) {
                //grow protein sequence buffer
                iProteinBuffer *= 2;
                if (0 == (sProteinSequence = (char*) realloc(sProteinSequence, iProteinBuffer*sizeof(char)))) {
                    fprintf(stderr, "\nERROR\tCould not allocate memory for protein %s (length %d)\n", sProteinReference, iProteinBuffer);
                    tempest_exit(EXIT_FAILURE);
                }
            }

            // concatenate to sequence
            strncat(sProteinSequence,sLineBuffer,iLengthBuffer);
            iLengthProtein += iLengthBuffer;
        }
    }

    // check for errors
    if (ferror(fp)) {
        fprintf(stderr, "\nERROR\tUnable to read from %s, while reading database.\n", args.sFasta);
        tempest_exit(EXIT_FAILURE);
    }
    
    // Final protein
    strncpy(&sProteinReferences[iProtein * MAX_LENGTH_REFERENCE],sProteinReference,MAX_LENGTH_REFERENCE);
    digest_protein(iProtein++, iLengthProtein, sProteinSequence);
    
    if (args.iPrintLevel && !args.bPrintCandidates) {
        fProgress = (100.0f * iProtein) / tempest.iNumProteins;
        printf("\b\b\b\b%*d%%", 3, (int) fProgress);
    }
        
    // Score Remaining Candidates
    for(i=0;i<tempest.iNumPrecursorBins;i++) {
        for (e = eScanIndex[i]; e; e = e->next) {
            //printf("%d\t%d\n", gpu_info.iNumScoringKernels, i);
            if (0 < e->iNumBufferedCandidates) gpu_score_candidates(e);
        }
    }
    
    // close fasta
    fclose(fp);

    // summary
    if (args.iPrintLevel && !args.bPrintCandidates) {
        printf("\n");
    }

    //cleanup
    free(sProteinSequence);
}

/*
 *  Digest protein sequence.
 */

void digest_protein(int iProtein, int iLengthProtein, char *sProteinSequence)
{
    static int i, iStart, iEnd, iPeptideLength, iEndLength, iNextNterm, iAA, iBin, iPattern;
    static int iNumAAModTargets, iModMass, iObjMass, iNumMods;
    static double dPeptideMass;
    static double dBaseMass;
    static float fModMass;
    static unsigned int uiMod;
    static unsigned int uiMaxModPattern;
    static bool bMatched;
    static eObj* e;
    static bool* bDigestSites = 0;
    static bool* bDigestNoSites = 0;
    static bool bNterm;
    static int  iDigestSites;
    static int modInds[MAX_MOD_TARGETS_PER_PEPTIDE];
    
    //setup enzyme sites
    if (bDigestSites == 0) {
        bDigestSites = (bool*) calloc(26, sizeof(bool));
        bDigestNoSites = (bool*) calloc(26, sizeof(bool));
        
        for (i=0; i<strlen(params.sDigestSites); i++) {
            bDigestSites[aa2i(params.sDigestSites[i])] = 1;
        }

        for (i=0; i<strlen(params.sDigestNoSites); i++) {
            bDigestNoSites[aa2i(params.sDigestNoSites[i])] = 1;
        }
    }

    // Initialize some variables.
    dBaseMass = params.dPeptideNtermMass + params.dPeptideCtermMass + H_MASS + OH_MASS + PROTON_MASS;
    iEndLength = iLengthProtein - 1;
    iStart = 0;
    bNterm = 1;

    // loop through protein
    while (iStart < iLengthProtein) {
        // reset
        iEnd = iStart;
        dPeptideMass = dBaseMass;
        if (iStart==0) dPeptideMass += params.dProteinNtermMass;
        iDigestSites = 0;
        iNextNterm = 0;
        iNumAAModTargets = 0;

        for (iEnd=iStart; iEnd < iLengthProtein; iEnd++) {
            iAA = aa2i(sProteinSequence[iEnd]);
            iPeptideLength = iEnd - iStart + 1;

            // APPLY DIGESTION RULES
            // What is the optimal arrangement for these checks?

            // too long?
            if (iPeptideLength > params.iMaxPeptideLength) {
                if (args.iPrintLevel==4) printf(">>> TOO LONG %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iPeptideLength, params.iMaxPeptideLength);
                break;
            }
            
            // update mass
            dPeptideMass += dMassAA[iAA];
            if (iEnd == iEndLength) dPeptideMass += params.dProteinCtermMass;
            
            // too large (no spectra to match it)?
            if (dPeptideMass > tempest.dMaxPeptideMass) {
                if (args.iPrintLevel==4) printf(">>> TOO LARGE %d %.*s (%.4f > %.4f)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, dPeptideMass, tempest.dMaxPeptideMass);
                break;
            }
            
            // mod target?
            if (cModSites[iAA]) {
                modInds[iNumAAModTargets] = iEnd - iStart;
                iNumAAModTargets++;
                
                // too many mod targets?
                if (iNumAAModTargets > MAX_MOD_TARGETS_PER_PEPTIDE) {
                    if (args.iPrintLevel==4) printf(">>> TOO MANY MOD TARGETS %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iNumAAModTargets, MAX_MOD_TARGETS_PER_PEPTIDE);
                    break;
                }
            }
            
            // cleavage site?
            if ((iEnd == iEndLength) || (bDigestSites[iAA] && !bDigestNoSites[aa2i(sProteinSequence[(iEnd+1)])])) {
                
                // set next iStart value
                if (!iNextNterm) {
                    iNextNterm = iEnd + 1;
                }

                // too many missed cleavages?
                if (iDigestSites > params.iDigestMaxMissedCleavages) {
                    if (args.iPrintLevel==4) printf(">>> TOO MANY MISSED CLEAVAGES %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iDigestSites, params.iDigestMaxMissedCleavages);
                    break;
                }

                iDigestSites++;
            }
            else if (params.iDigestSpecificity == 2 || (params.iDigestSpecificity == 1 && bNterm == 0)) {
                if (args.iPrintLevel==4) printf(">>> NONSPECIFIC %d %.*s (%d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, params.iDigestSpecificity);
                continue;
            }
            
            // too short or too small?
            if (iPeptideLength < params.iMinPeptideLength) {
                if (args.iPrintLevel==4) printf(">>> TOO SHORT %d %.*s (%d < %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iPeptideLength, params.iMinPeptideLength);
                continue;
            }

            //too small?
            if (dPeptideMass < tempest.dMinPeptideMass) {
                if (args.iPrintLevel==4) printf(">>> TOO SMALL %d %.*s (%.4f < %.4f)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, dPeptideMass, tempest.dMinPeptideMass);
                continue;
            }
            
            // Potential Peptide found
            tempest.lNumPeptides += 1;

            if (args.bNoMatch) {
                continue;
            }
            
            // Loop spectra for matches
            iObjMass = (int) roundf((float) dPeptideMass / tempest.fPrecursorBinWidth);
            uiMaxModPattern = 1 << iNumAAModTargets;

            cObj init_cCandidate;
            strncpy(init_cCandidate.sPeptide, &sProteinSequence[iStart], (size_t)iPeptideLength);
            init_cCandidate.sPeptide[iPeptideLength] = '\0';
            init_cCandidate.bNtermMod = 0;
            init_cCandidate.bCtermMod = 0;
            init_cCandidate.iProtein       = iProtein;
            init_cCandidate.iPeptideLength = iPeptideLength;
            init_cCandidate.fPeptideMass   = float(dPeptideMass);
            init_cCandidate.cBefore        = (iStart == 0)        ? '-' : sProteinSequence[iStart-1];
            init_cCandidate.cAfter         = (iEnd == iEndLength) ? '-' : sProteinSequence[iEnd+1];

            gen_candidates(init_cCandidate, params.cVariableNtermSymbol, params.cVariableCtermSymbol, modInds, iNumAAModTargets, 0);
        }

        //find next n-terminal residue

        // Full: next digest site
        if (params.iDigestSpecificity == 2) {
            if (iNextNterm > 0) {
                iStart = iNextNterm;
            }
            else {
                for (iStart=iEnd+1;
                     iStart < iLengthProtein && (!bDigestSites[aa2i(sProteinSequence[(iStart-1)])] || bDigestNoSites[aa2i(sProteinSequence[iStart])]);
                     iStart++);
            }
        }

        // Partial or Nonspecific: next residue
        else {
            iStart++;
            if (bDigestSites[aa2i(sProteinSequence[(iStart-1)])] && !bDigestNoSites[aa2i(sProteinSequence[iStart])]) {
                bNterm = 1;
            }
        }
    }
}
    
/*
 * Count reference lines in a fasta file.
 */
int count_references(FILE *fp) {
    int n = 0;
    char c;
    char cPrevious = '\n';
    
    while (EOF != (c = getc (fp))) {
        if (c == '>' && cPrevious == '\n') n++;
        cPrevious = c;
    }
    rewind(fp);
    return n;
}

/*
 * Parse a (truncated) protein reference from a fasta description line.
 *
 *  e.g.
 *  description line:  >P31946|1433B_HUMAN 14-3-3 protein beta/alpha - Homo sapiens
 *  protein reference: P31946|1433B_HUMAN
 *
 *
 *      sLine       null-terminated protein description from fasta file (with a leading '>').
 *      sReference  pointer to pre-allocated memory for the parsed protein reference.
 */

void parse_reference(char* sLine, char* sReference) {
    static int iLengthRef;
    
    //look for first space
    for (iLengthRef=1; iLengthRef<MAX_LENGTH_REFERENCE && sLine[iLengthRef+1]!=' '; iLengthRef++);

    //copy
    strncpy(sReference, sLine+1, iLengthRef);
    sReference[iLengthRef] = '\0';
}

void gen_candidates(cObj cCandidate, bool ntermMod, bool ctermMod, int* modInds, int modIndsLeft, int modCount) {
    if (modCount >= params.iModificationsMax)
        store_candidate(cCandidate);
    else if (ntermMod) {
        cObj new_cCandidate;
        memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
        new_cCandidate.bNtermMod = 1;
        new_cCandidate.fPeptideMass += params.dVariableNtermMassDiff;
        gen_candidates(cCandidate, 0, ctermMod, modInds, modIndsLeft, modCount);
        gen_candidates(new_cCandidate, 0, ctermMod, modInds, modIndsLeft, modCount+1);
    }
    else if (modIndsLeft > 0) {
        int modInd = *modInds;
        int aai = aa2i(cCandidate.sPeptide[modInd]);
        cObj new_cCandidate;
        memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
        new_cCandidate.sPeptide[modInd] = tolower(new_cCandidate.sPeptide[modInd]);
        new_cCandidate.fPeptideMass += fModValues[aai];
        gen_candidates(cCandidate, ntermMod, ctermMod, modInds+1, modIndsLeft-1, modCount);
        gen_candidates(new_cCandidate, ntermMod, ctermMod, modInds+1, modIndsLeft-1, modCount+1);
    }
    else if (ctermMod) {
        cObj new_cCandidate;
        memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
        new_cCandidate.bCtermMod = 1;
        new_cCandidate.fPeptideMass += params.dVariableCtermMassDiff;
        gen_candidates(cCandidate, ntermMod, 0, modInds, modIndsLeft, modCount);
        gen_candidates(new_cCandidate, ntermMod, 0, modInds, modIndsLeft, modCount+1);
    }
    else 
        store_candidate(cCandidate);
}

void store_candidate(cObj cCandidate) {
    float fModMass = cCandidate.fPeptideMass;
    int iModMass = (int) roundf(fModMass / tempest.fPrecursorBinWidth);
    for (int iBin=iModMass-1; iBin<=iModMass+1 && iBin<=tempest.iNumPrecursorBins; iBin++) {
        for (eObj* e = eScanIndex[iBin]; e != 0; e = e->next) {
            // check mass error
            if (params.bPrecursorTolerancePPM) {
                if (1000000.0f*fabs(fModMass - (float) e->dPrecursorMass) / fModMass > params.fPrecursorTolerance)
                    continue;
            }
            else {
                if (fabs(fModMass - (float) e->dPrecursorMass) > params.fPrecursorTolerance)
                    continue;
            }
                  
            tempest.lNumPSMs += 1;
            if (e->iNumBufferedCandidates == 0) {
                clWaitForEvents(1, &(e->clEventSent));
                if (PROFILE) {
                    cl_ulong start;
                    cl_ulong end;
                    int err;
                    err = clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
                    err |= clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
                    if (err == 0)
                        totalSendTime += (end-start);
                }
                clReleaseEvent(e->clEventSent);
            }
            e->pCandidateBufferFill[e->iNumBufferedCandidates] = cCandidate;
            e->iNumCandidates++;
            e->iNumBufferedCandidates++;
            if (e->iNumBufferedCandidates == config.iCandidateBufferSize) {
                //printf("%d\t%d\n", gpu_info.iNumScoringKernels, iBin);
                gpu_score_candidates(e);
            }
        }
    }
}
