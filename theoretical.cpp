#include "tempest.hpp"

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

extern void Tempest::search_fasta_database()
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
    float fProgress;

    if (Tempest::args.bNoDigest) return;

    // Setup

    // allocate working memory
    if (0 == (sProteinSequence = (char*) malloc(iProteinBuffer*sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for the active protein sequence.\n");
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // calculate peptide mass range
    float fMaxAAModMass = 0.0f;
    for (i=0; i<Tempest::params.iNumMods; i++) {
        if (Tempest::params.tMods[i].cSymbol && Tempest::params.tMods[i].dMassDiff > fMaxAAModMass)
            fMaxAAModMass = Tempest::params.tMods[i].dMassDiff;
    }
    float fMaxNtermModMass = 0.0f;
    for (i=0; i<Tempest::numNtermModSites; i++)
        if (ntermModMasses[i] > fMaxNtermModMass)
            fMaxNtermModMass = ntermModMasses[i];
    float fMaxCtermModMass = 0.0f;
    for (i=0; i<Tempest::numCtermModSites; i++)
        if (Tempest::ctermModMasses[i] > fMaxCtermModMass)
            fMaxCtermModMass = Tempest::ctermModMasses[i];

    float maxMaxTermModMass = fMaxNtermModMass > fMaxCtermModMass ? fMaxNtermModMass : fMaxCtermModMass;
    float minMaxTermModMass = fMaxNtermModMass < fMaxCtermModMass ? fMaxNtermModMass : fMaxCtermModMass;

    float fMaxMassDiff = 0.0;
    int remainingMods = Tempest::params.iModificationsMax;
    if (remainingMods && maxMaxTermModMass > fMaxAAModMass) {
        fMaxMassDiff += maxMaxTermModMass;
        remainingMods -= 1;
    }
    if (remainingMods && minMaxTermModMass > fMaxAAModMass) {
        fMaxMassDiff += minMaxTermModMass;
        remainingMods -= 1;
    }
    fMaxMassDiff += fMaxAAModMass * remainingMods;

    if (Tempest::params.bPrecursorTolerancePPM) {
        Tempest::tempest.fMinPrecursorMass -= Tempest::tempest.fMinPrecursorMass*Tempest::params.fPrecursorTolerance/1000000.0f;
        Tempest::tempest.fMaxPrecursorMass += Tempest::tempest.fMaxPrecursorMass*Tempest::params.fPrecursorTolerance/1000000.0f;
    }
    else {
        Tempest::tempest.fMinPrecursorMass -= Tempest::params.fPrecursorTolerance;
        Tempest::tempest.fMaxPrecursorMass += Tempest::params.fPrecursorTolerance;
    }
    Tempest::tempest.dMinPeptideMass = Tempest::tempest.fMinPrecursorMass - fMaxMassDiff;
    Tempest::tempest.dMaxPeptideMass = Tempest::tempest.fMaxPrecursorMass;
    
    // open fasta file
    if (0 == (fp = (FILE *) fopen(Tempest::args.sFasta, "r"))) {
        fprintf(stderr, "\nERROR\tUnable to open fasta database: %s\n", Tempest::args.sFasta);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // count proteins
    Tempest::tempest.iNumProteins = count_references(fp);
    
    if (Tempest::tempest.iNumProteins == 0) {
        fprintf(stderr, "\nERROR\tNo proteins found in fasta database: %s\n", Tempest::args.sFasta);
        Tempest::tempest_exit(EXIT_FAILURE);
    }

    // progress
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
        printf(" » Digesting %d proteins...     ", Tempest::tempest.iNumProteins);
    }
    
    // Get memory for the protein references
    if (0 == (Tempest::sProteinReferences = (char *) malloc(Tempest::tempest.iNumProteins * MAX_LENGTH_REFERENCE * sizeof(char)))) {
        fprintf(stderr, "\nERROR\tUnable to allocate memory for protein references.\n");
        Tempest::tempest_exit(EXIT_FAILURE);
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
                strcpy(&Tempest::sProteinReferences[iProtein * MAX_LENGTH_REFERENCE * sizeof(char)], sProteinReference);
                digest_protein(iProtein++, iLengthProtein, sProteinSequence);

                if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
                    fProgress = (100.0f * iProtein) / Tempest::tempest.iNumProteins;
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
                    Tempest::tempest_exit(EXIT_FAILURE);
                }
            }

            // concatenate to sequence
            strncat(sProteinSequence,sLineBuffer,iLengthBuffer);
            iLengthProtein += iLengthBuffer;
        }
    }

    // check for errors
    if (ferror(fp)) {
        fprintf(stderr, "\nERROR\tUnable to read from %s, while reading database.\n", Tempest::args.sFasta);
        Tempest::tempest_exit(EXIT_FAILURE);
    }
    
    // Final protein
    strncpy(&Tempest::sProteinReferences[iProtein * MAX_LENGTH_REFERENCE],sProteinReference,MAX_LENGTH_REFERENCE);
    digest_protein(iProtein++, iLengthProtein, sProteinSequence);
    
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
        fProgress = (100.0f * iProtein) / Tempest::tempest.iNumProteins;
        printf("\b\b\b\b%*d%%", 3, (int) fProgress);
    }
        
    // Score Remaining Candidates
    for(i=0;i<Tempest::tempest.iNumPrecursorBins;i++) {
        for (e = Tempest::eScanIndex[i]; e; e = e->next) {
            //printf("%d\t%d\n", gpu_info.iNumScoringKernels, i);
            if (e->iNumBufferedCandidates > 0)
                e->device->scoreCandidates(e);
        }
    }
    
    // close fasta
    fclose(fp);

    // summary
    if (Tempest::args.iPrintLevel && !Tempest::args.bPrintCandidates) {
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
    int iStart, iEnd, iPeptideLength, iEndLength, iNextNterm, iBin, iPattern;
    int iNumAAModTargets, iModMass, iNumMods;
    double dPeptideMass;
    double dBaseMass;
    unsigned int uiMod;
    unsigned int uiMaxModPattern;
    bool bMatched;
    eObj* e;
    bool bNterm;
    int  iDigestSites;
    int modInds[MAX_MOD_TARGETS_PER_PEPTIDE];

    // Initialize some variables.
    dBaseMass = Tempest::params.dPeptideNtermMass + Tempest::params.dPeptideCtermMass + H_MASS + OH_MASS + PROTON_MASS;
    iEndLength = iLengthProtein - 1;
    iStart = 0;
    bNterm = 1;

    // loop through protein
    while (iStart < iLengthProtein) {
        // reset
        iEnd = iStart;
        dPeptideMass = dBaseMass;
        if (iStart==0) dPeptideMass += Tempest::params.dProteinNtermMass;
        iDigestSites = 0;
        iNextNterm = 0;
        iNumAAModTargets = 0;

        for (iEnd=iStart; iEnd < iLengthProtein; iEnd++) {
            unsigned char aa = sProteinSequence[iEnd];
            
            iPeptideLength = iEnd - iStart + 1;

            // APPLY DIGESTION RULES
            // What is the optimal arrangement for these checks?

            // too long?
            if (iPeptideLength >= Tempest::params.iMaxPeptideLength) {
                // if (Tempest::args.iPrintLevel==4) printf(">>> TOO LONG %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iPeptideLength, Tempest::params.iMaxPeptideLength);
                break;
            }
            
            // update mass
            dPeptideMass += Tempest::dMassAA[aa];
            if (iEnd == iEndLength) dPeptideMass += Tempest::params.dProteinCtermMass;
            
            // too large (no spectra to match it)?
            if (dPeptideMass > Tempest::tempest.dMaxPeptideMass) {
                // if (Tempest::args.iPrintLevel==4) printf(">>> TOO LARGE %d %.*s (%.4f > %.4f)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, dPeptideMass, Tempest::tempest.dMaxPeptideMass);
                break;
            }
            
            // mod target?
            if (Tempest::numAAModSites[aa]) {
                // too many mod targets?
                if (iNumAAModTargets >= MAX_MOD_TARGETS_PER_PEPTIDE) {
                    // if (Tempest::args.iPrintLevel==4)
                    //     printf(">>> TOO MANY MOD TARGETS %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iNumAAModTargets, MAX_MOD_TARGETS_PER_PEPTIDE);
                    break;
                }
                modInds[iNumAAModTargets] = iEnd - iStart;
                iNumAAModTargets++;
            }
            
            // cleavage site?
            if ((iEnd == iEndLength) || (Tempest::bDigestSites[aa] && !Tempest::bDigestNoSites[sProteinSequence[(iEnd+1)]])) {
                
                // set next iStart value
                if (!iNextNterm) {
                    iNextNterm = iEnd + 1;
                }

                // too many missed cleavages?
                if (iDigestSites > Tempest::params.iDigestMaxMissedCleavages) {
                    // if (Tempest::args.iPrintLevel==4) printf(">>> TOO MANY MISSED CLEAVAGES %d %.*s (%d > %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iDigestSites, Tempest::params.iDigestMaxMissedCleavages);
                    break;
                }

                iDigestSites++;
            }
            else if (Tempest::params.iDigestSpecificity == 2 || (Tempest::params.iDigestSpecificity == 1 && bNterm == 0)) {
                // if (Tempest::args.iPrintLevel==4) printf(">>> NONSPECIFIC %d %.*s (%d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, Tempest::params.iDigestSpecificity);
                continue;
            }
            
            // too short or too small?
            if (iPeptideLength < Tempest::params.iMinPeptideLength) {
                // if (Tempest::args.iPrintLevel==4) printf(">>> TOO SHORT %d %.*s (%d < %d)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, iPeptideLength, Tempest::params.iMinPeptideLength);
                continue;
            }

            //too small?
            if (dPeptideMass < Tempest::tempest.dMinPeptideMass) {
                // if (Tempest::args.iPrintLevel==4) printf(">>> TOO SMALL %d %.*s (%.4f < %.4f)\n", iProtein, iEnd-iStart+1, sProteinSequence+iStart, dPeptideMass, Tempest::tempest.dMinPeptideMass);
                continue;
            }
            
            // Potential Peptide found
            Tempest::tempest.lNumPeptides += 1;

            // if (Tempest::args.bNoMatch) {
            //     continue;
            // }

            cObj init_cCandidate;
            strncpy((char*)init_cCandidate.sPeptide, &sProteinSequence[iStart], (size_t)iPeptideLength);
            init_cCandidate.sPeptide[iPeptideLength] = '\0';
            init_cCandidate.ntermMod = 0;
            init_cCandidate.ctermMod = 0;
            init_cCandidate.iProtein       = iProtein;
            init_cCandidate.iPeptideLength = iPeptideLength;
            init_cCandidate.fPeptideMass   = float(dPeptideMass);
            init_cCandidate.cBefore        = (iStart == 0)        ? '-' : sProteinSequence[iStart-1];
            init_cCandidate.cAfter         = (iEnd == iEndLength) ? '-' : sProteinSequence[iEnd+1];

            gen_candidates(init_cCandidate, Tempest::numNtermModSites>0, Tempest::numCtermModSites>0, modInds, iNumAAModTargets, 0);
        }

        //find next n-terminal residue

        // Full: next digest site
        if (Tempest::params.iDigestSpecificity == 2) {
            if (iNextNterm > 0) {
                iStart = iNextNterm;
            }
            else {
                for (iStart=iEnd+1;
                     iStart < iLengthProtein && (!Tempest::bDigestSites[sProteinSequence[(iStart-1)]] || Tempest::bDigestNoSites[sProteinSequence[iStart]]);
                     iStart++);
            }
        }

        // Partial or Nonspecific: next residue
        else {
            iStart++;
            bNterm = Tempest::bDigestSites[sProteinSequence[(iStart-1)]] && !Tempest::bDigestNoSites[sProteinSequence[iStart]];
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
    if (modCount >= Tempest::params.iModificationsMax)
        store_candidate(cCandidate);
    else if (modIndsLeft > 0) {
        int modInd = *modInds;
        unsigned char aa = cCandidate.sPeptide[modInd];
        for (int i=0; i<Tempest::numAAModSites[aa]; i++) {
            cObj new_cCandidate = cCandidate;
            //memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            new_cCandidate.sPeptide[modInd] = Tempest::toMod(aa, i+1);
            new_cCandidate.fPeptideMass += Tempest::fModValues[aa][i];
            //printf("%c[%d] = %f\n", aa, i, sPeptide[modInd], Tempest::fModValues[aa][i]);
            gen_candidates(new_cCandidate, ntermMod, ctermMod, modInds+1, modIndsLeft-1, modCount+1);
        }
        gen_candidates(cCandidate, ntermMod, ctermMod, modInds+1, modIndsLeft-1, modCount);
    }
    else if (ntermMod) {
        cObj new_cCandidate;
        for (int i=0; i<Tempest::numNtermModSites; i++) {
            memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            new_cCandidate.ntermMod = i+1;
            new_cCandidate.fPeptideMass += Tempest::ntermModMasses[i];
            gen_candidates(new_cCandidate, 0, ctermMod, modInds, modIndsLeft, modCount+1);
        }
        gen_candidates(cCandidate, 0, ctermMod, modInds, modIndsLeft, modCount);
    }
    else if (ctermMod) {
        cObj new_cCandidate;
        for (int i=0; i<Tempest::numCtermModSites; i++) {
            memcpy(&new_cCandidate, &cCandidate, sizeof(cObj));
            new_cCandidate.ctermMod = i+1;
            new_cCandidate.fPeptideMass += Tempest::ctermModMasses[i];  
            gen_candidates(new_cCandidate, ntermMod, 0, modInds, modIndsLeft, modCount+1);
        }
        gen_candidates(cCandidate, ntermMod, 0, modInds, modIndsLeft, modCount);
    }
    else 
        store_candidate(cCandidate);
}

void store_candidate(cObj cCandidate) {
    float fModMass = cCandidate.fPeptideMass;
    if (fModMass < Tempest::tempest.fMinPrecursorMass || fModMass > Tempest::tempest.fMaxPrecursorMass)
        return;
    int iModMass = (int) roundf(fModMass / Tempest::tempest.fPrecursorBinWidth);
    for (int iBin=iModMass-1; iBin<=iModMass+1 && iBin<=Tempest::tempest.iNumPrecursorBins; iBin++) {
        for (eObj* e = Tempest::eScanIndex[iBin]; e != 0; e = e->next) {
            // check mass error
            if (Tempest::params.bPrecursorTolerancePPM) {
                if (1000000.0f*fabs(fModMass - (float) e->dPrecursorMass) / fModMass > Tempest::params.fPrecursorTolerance)
                    continue;
            }
            else {
                if (fabs(fModMass - (float) e->dPrecursorMass) > Tempest::params.fPrecursorTolerance)
                    continue;
            }
                  
            Tempest::tempest.lNumPSMs += 1;
            if (e->iNumBufferedCandidates == 0) {
                clWaitForEvents(1, &(e->clEventSent));
                /*
                if (PROFILE) {
                    cl_ulong start;
                    cl_ulong end;
                    int err;
                    err = clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &start, NULL);
                    err |= clGetEventProfilingInfo(e->clEventSent, CL_PROFILING_COMMAND_END,   sizeof(cl_ulong), &end,   NULL);
                    if (err == 0)
                        totalSendTime += (end-start);
                }
                */
                clReleaseEvent(e->clEventSent);
            }
            e->candidateBuffer[e->iNumBufferedCandidates] = cCandidate;
            //memcpy(e->candidateBuffer+e->iNumBufferedCandidates, &cCandidate, sizeof(cObj));
            e->iNumCandidates++;
            e->iNumBufferedCandidates++;
            if (e->iNumBufferedCandidates == e->candidateBufferSize) {
                //printf("%d\t%d\n", gpu_info.iNumScoringKernels, iBin);
                e->device->scoreCandidates(e);
            }
        }
    }
}
